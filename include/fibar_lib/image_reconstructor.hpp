// -*-c++-*---------------------------------------------------------------------------------------
// Copyright 2025 Bernd Pfrommer <bernd.pfrommer@gmail.com>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef FIBAR_LIB_IMAGE_RECONSTRUCTOR_HPP
#define FIBAR_LIB_IMAGE_RECONSTRUCTOR_HPP

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <deque>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#include "fibar_lib/spatial_filter.hpp"
#include "fibar_lib/state.hpp"

#define LIKELY(x) __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)
#define SANITY_CHECKS

namespace fibar_lib
{
template <int tile_size>
inline size_t getTileIndex(uint16_t ex, uint16_t ey, uint16_t tile_stride_y)
{
  return ((ey / tile_size) * tile_stride_y + (ex / tile_size) * tile_size);
}
template <>
inline size_t getTileIndex<2>(uint16_t ex, uint16_t ey, uint16_t tile_stride_y)

{
  return ((ey >> 1) * tile_stride_y + (ex & ~1));
}

template <bool filter_spatially>
class BaseImageReconstructor
{
public:
  using state_t = float;

  void initialize(size_t width, size_t height, uint32_t cutoff_time)
  {
    width_ = width;
    height_ = height;
    // compute filter coefficients
    double alpha(0);
    double beta(0);
    computeAlphaBeta(static_cast<double>(cutoff_time), &alpha, &beta);
    c_[0] = static_cast<float>(alpha);
    c_[1] = static_cast<float>(1.0 - alpha);
    c_[2] = static_cast<float>(beta);
    c_[3] = static_cast<float>(0.5 * (1 + beta));
    state_.resize(width * height, State<filter_spatially>());
  }

  inline void update_filter(
    State<filter_spatially> & s, uint32_t t, uint16_t ex, uint16_t ey, uint8_t polarity)
  {
    const auto p = static_cast<state_t>((polarity == 0) ? -1 : 1);
#ifdef RESCALE
    // change in polarity, will be scale * (0 or +-2)
    const auto dp = s.getScale() * static_cast<float>(p - s.getPbar());
#else
    // raw change in polarity, will be 0 or +-2
    const auto dp = static_cast<float>(p - s.getPbar());
#endif
    // run the temporal filter
    const auto L = c_[2] * s.getL() + c_[3] * dp;
    // update state
    s.setPbar(s.getPbar() * c_[0] + p * c_[1]);
    s.setL(L);
  }
  size_t getWidth() const { return (width_); }
  size_t getHeight() const { return (height_); }
  const auto & getState() const { return (state_); }

  void computeAlphaBeta(const double T_cut, double * alpha, double * beta)
  {
    // compute the filter coefficients alpha and beta (see frequency cam paper)
    const double omega_cut = 2 * M_PI / T_cut;
    const double phi = 2 - std::cos(omega_cut);
    *alpha = (1.0 - std::sin(omega_cut)) / std::cos(omega_cut);
    *beta = phi - std::sqrt(phi * phi - 1.0);  // see paper
  }

  void getImage(uint8_t * img, size_t stride) const
  {
    // find min and max for normalization
    float min_L = std::numeric_limits<float>::max();
    float max_L = std::numeric_limits<float>::min();
    for (size_t i = 0; i < height_ * width_; i++) {
      if (state_[i].getL() > max_L) {
        max_L = state_[i].getL();
      }
      if (state_[i].getL() < min_L) {
        min_L = state_[i].getL();
      }
    }
    // copy image over

    const float scale = 255.0F / (max_L - min_L);
    for (size_t iy = 0; iy < height_; iy++) {
      const size_t y_off = iy * stride;
      const size_t y_off_state = iy * width_;
      for (size_t ix = 0; ix < width_; ix++) {
        const auto & s = state_[y_off_state + ix];
        img[y_off + ix] = static_cast<uint8_t>((s.getL() - min_L) * scale);
      }
    }
  }

  // ------------------- variables ------------------
  size_t width_{0};
  size_t height_{0};
  std::vector<State<filter_spatially>> state_;  // filter state
  std::array<float, 4> c_{0, 0, 0, 0};          // filter coefficients
};

// ---------------------------------------------------------------------
// template declaration
//
template <bool filter_spatially = true, uint8_t tile_size = 2>
class ImageReconstructor : public BaseImageReconstructor<filter_spatially>
{
public:
  ImageReconstructor() = default;
};

// ---------------------------------------------------------------------
// specialization for no spatial filtering
//
template <uint8_t tile_size>
class ImageReconstructor<false, tile_size> : public BaseImageReconstructor<false>
{
public:
  ImageReconstructor() = default;
  void event(uint32_t t, uint16_t ex, uint16_t ey, uint8_t polarity)
  {
    auto & s = state_[ey * width_ + ex];
    BaseImageReconstructor<false>::update_filter(s, t, ex, ey, polarity);
  }
};
// ---------------------------------------------------------------------
// specialization for spatial filtering
//
template <uint8_t tile_size>
class ImageReconstructor<true, tile_size> : public BaseImageReconstructor<true>
{
public:
  ImageReconstructor() = default;
  using state_t = float;
  static constexpr std::array<std::array<state_t, 3>, 3> GAUSSIAN_3x3 = {
    {{0.0625, 0.125, 0.0625}, {0.125, 0.25, 0.125}, {0.0625, 0.125, 0.0625}}};
  static constexpr std::array<std::array<state_t, 5>, 5> GAUSSIAN_5x5 = {
    {{0.003663, 0.01465201, 0.02564103, 0.01465201, 0.003663},
     {0.01465201, 0.05860806, 0.0952381, 0.05860806, 0.01465201},
     {0.02564103, 0.0952381, 0.15018315, 0.0952381, 0.02564103},
     {0.01465201, 0.05860806, 0.0952381, 0.05860806, 0.01465201},
     {0.003663, 0.01465201, 0.02564103, 0.01465201, 0.003663}}};

  void event(uint32_t t, uint16_t ex, uint16_t ey, uint8_t polarity)
  {
    auto & s = state_[ey * width_ + ex];
    BaseImageReconstructor<true>::update_filter(s, t, ex, ey, polarity);
    if (!s.isActive()) {
      num_occupied_pixels_ += fill_ratio_denom_;
      // state of top left corner of tile has actual pixel-in-tile count
      auto & tile = state_[getTileIdx(ex, ey)];
      if (tile.getNumPixActive() == 0) {
        num_occupied_tiles_ += fill_ratio_num_;  // first active pixel in this tile
      }
      tile.incNumPixActive();  // bump number of pixels in this tile
    }
    s.incNumEventsInQueue();
    events_.push_back(Event(ex, ey, static_cast<int8_t>(polarity)));
    processEventQueue();  // adjusts size of event window
  };

  size_t getCurrentQueueSize() const { return (events_.size()); }
  double getCurrentFillRatio() const
  {
    return (
      num_occupied_tiles_ == 0 ? -1.0
                               : static_cast<double>(num_occupied_pixels_) /
                                   (num_occupied_tiles_ * tile_size * tile_size));
  }

  inline size_t getEventWindowSize() const { return (event_window_size_); }
  void initialize(size_t width, size_t height, uint32_t cutoff_time, double fill_ratio)
  {
    BaseImageReconstructor<true>::initialize(width, height, cutoff_time);
    tile_stride_y_ = width * tile_size;
    constexpr int max_area = 1 << 7;
    if (tile_size * tile_size > max_area) {
      // guard against overflow of count of occupied pixels in tile
      std::cerr << "activity tile size too big: " << tile_size << " must be < "
                << static_cast<int>(std::sqrt(max_area)) << std::endl;
      throw(std::runtime_error("activity tile size too big"));
    }
    // disable any queue usage if tile size is set to zero
    max_window_size_ = tile_size > 0 ? static_cast<uint64_t>(width_ * height_) : 0;
    setFillRatio(fill_ratio);
  }

  void getActivePixelImage(uint8_t * img, size_t stride) const
  {
    // clear image
    memset(img, 0, height_ * stride);
    for (const auto & qe : events_) {
      img[qe.y() * stride + qe.x()]++;
    }
  }

  void setFillRatio(double fill_ratio)
  {
    fill_ratio_denom_ = 100;
    // A is the area of the tile (in pixels)
    const double A = static_cast<double>(tile_size * tile_size);
    // how many tiles per pixel when fully filled
    const double tiles_per_pixel = 1.0 / A;
    // a fill ratio below 1 pixel per tile is not achievable
    const double r = std::min(1.0, std::max(fill_ratio, tiles_per_pixel + 1e-3));
    const double np_nt = A * r;  // targeted number of pixels per tile
    fill_ratio_num_ = static_cast<uint64_t>(np_nt * fill_ratio_denom_);
    // The update equation for the queue length q is
    // q_{k+1} = floor(q_k * f)
    // where f is the current gain:
    // f = (num_tiles * 100 * np_nt) / (num_pixels * 100)
    // For the queue to be able to grow, we must ensure that
    // q_k * f > q_k + 1
    // meaning q_k > 1/(f - 1)
    // The largest that f can become is when num_tiles == num_pixels,
    // in which case f = np_t, and so q_k > 1 / (np_nt - 1)
    //
    min_window_size_ = A > 0 ? std::ceil((1.0 / (np_nt - 1.0))) : 0;
  }

#ifdef RESCALE
  void readScaleFile(const std::string & fname)
  {
    std::ifstream file;
    file.open(fname);
    if (!file.is_open()) {
      throw std::runtime_error("cannot open scale file: " + fname);
    }
    const uint32_t n_pix = width_ * height_;
    uint64_t sum{0};
    uint32_t s;
    for (size_t idx = 0; (file >> s) && (idx < n_pix * 2); idx++) {
      sum += s;
    }
    const double ntot_avg = sum / n_pix;
    file.close();
    file.open(fname);
    uint32_t n_on, n_off;
    size_t idx = 0;
    double ss{0}, ss2{0}, sum_inv{0};
    for (; (file >> n_on) && (file >> n_off) && (idx < n_pix); idx++) {
      const double C_i = ntot_avg / static_cast<double>(n_on + n_off);
      state_[idx].scale = C_i;
      ss += C_i;
      ss2 += C_i * C_i;
      sum_inv += 1.0 / C_i;
    }
    ss = ss / n_pix;
    ss2 = ss2 / n_pix;
    const double stddev = std::sqrt(ss2 - ss * ss);
    std::cout << "read scale file: " << fname << " with " << idx << " entries and " << ntot_avg
              << " events/pixel, avg C: " << ss << " stdev: " << stddev
              << " harmonic mean: " << (n_pix / sum_inv) << std::endl;
  }
#endif

private:
  inline size_t getTileIdx(uint16_t ex, uint16_t ey) const
  {
    return (getTileIndex<tile_size>(ex, ey, tile_stride_y_));
  }

  void processEventQueue()
  {
    while (events_.size() > event_window_size_) {
      const Event & e = events_.front();
      auto & s = state_[e.y() * width_ + e.x()];
#ifdef SANITY_CHECKS
      if (!s.isActive()) {
        std::cerr << "FIBAR: pixel at (" << e.x() << "," << e.y() << ") has bad activity counter!"
                  << std::endl;
        std::cerr << "FIBAR: likely this is a hot(bad) pixel. Mask it out in the camera driver!"
                  << std::endl;
        throw std::runtime_error(
          "bad activity counter at pixel " + std::to_string(e.x()) + "," + std::to_string(e.y()));
      }
#endif
      s.decNumEventsInQueue();
      if (!s.isActive()) {
#ifdef SPATIAL_FILTER_5x5
        s =
          spatial_filter::filter<State, 5>(&state_[0], e.x(), e.y(), width_, height_, GAUSSIAN_5x5);
#else
        // s =  spatial_filter::filter<State, 3>(&state_[0], e.x(), e.y(), width_, height_, GAUSSIAN_3x3);
        s = spatial_filter::filter_3x3(state_.data(), e.x(), e.y(), width_, height_, GAUSSIAN_3x3);
#endif
        auto & tile = state_[getTileIdx(e.x(), e.y())];  // state of top left corner of tile
#ifdef SANITY_CHECKS
        if (tile.getNumPixActive() == 0) {
          std::cerr << e.x() << " " << e.y() << " tile " << getTileIdx(e.x(), e.y()) << " is empty!"
                    << std::endl;
          throw std::runtime_error(
            "empty tile at " + std::to_string(e.x()) + "," + std::to_string(e.y()));
        }
#endif
        // remove number of pixels in this tile
        tile.decNumPixActive();
        if (tile.getNumPixActive() == 0) {
          num_occupied_tiles_ -= fill_ratio_num_;
        }
        num_occupied_pixels_ -= fill_ratio_denom_;
      }
      events_.pop_front();  // remove element now
    }
    // adjust event window size up or down to match the fill ratio:
    // new_size = old_size * current_fill_ratio / desired_fill_ratio
    // The idea is that as the event window increases, the features will "fill out"
#define AVOID_DIVISION
#ifdef AVOID_DIVISION
    const int64_t ntfn = num_occupied_tiles_;
    int64_t npfd = num_occupied_pixels_;
    if (UNLIKELY(npfd <= 1)) {
      npfd = fill_ratio_denom_;
    }
    if (LIKELY(std::abs(500 * (ntfn - npfd)) > npfd)) {
      const uint64_t target_size = (event_window_size_ * ntfn) / npfd;
      // prevent the event window from collapsing to zero and from growing without bounds
      event_window_size_ = std::max(min_window_size_, std::min(max_window_size_, target_size));
    }
#else
    const uint64_t target_size = (event_window_size_ * num_occupied_tiles_) /
                                 (std::max(num_occupied_pixels_, fill_ratio_denom_));
    event_window_size_ = std::max(min_window_size_, std::min(max_window_size_, target_size));
#endif
  }

  class Event
  {
  public:
    explicit Event(uint16_t x = 0, uint16_t y = 0, int8_t p = 0) : ex(x), ey(y | (p << 15)) {}
    inline uint16_t x() const { return (ex); }
    inline uint16_t y() const { return (ey & 0x7fff); }
    inline int8_t p() const { return ((ey & 0x8000) >> 15); }

  private:
    uint16_t ex{0};
    uint16_t ey{0};
  };
  // ---------- related to activity detection
  static constexpr int START_WINDOW_SIZE = 2000;
  uint16_t tile_stride_y_{0};                      // size of stride in tiled image
  uint64_t event_window_size_{START_WINDOW_SIZE};  // current event window size
  uint64_t fill_ratio_denom_{2};                   // denominator of fill ratio
  uint64_t fill_ratio_num_{1};                     // numerator of fill ratio
  uint64_t num_occupied_pixels_{0};                // currently occupied number of pixels
  uint64_t num_occupied_tiles_{0};                 // currently occupied number of blocks
  uint64_t max_window_size_{0};                    // maximum size of event window
  uint64_t min_window_size_{0};                    // minimum size of event window
  std::deque<Event> events_;                       // queue with buffered events
};
}  // namespace fibar_lib
#endif  // FIBAR_LIB_IMAGE_RECONSTRUCTOR_HPP
