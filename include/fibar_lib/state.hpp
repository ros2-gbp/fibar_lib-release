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

#ifndef FIBAR_LIB_STATE_HPP
#define FIBAR_LIB_STATE_HPP

#include <array>
#include <cstddef>
#include <cstdint>

namespace fibar_lib
{
class BaseState
{
public:
  using state_t = float;
  explicit BaseState(state_t L_a = 0, state_t pbar_a = 0) : L(L_a), pbar(pbar_a) {}
  inline void operator+=(const BaseState & s)
  {
    L += s.L;
    // leave other fields untouched
  }

  inline BaseState operator*(const float c) const { return (BaseState(c * L)); }

  inline state_t getL() const { return (L); }
  inline state_t getPbar() const { return (pbar); }
  inline void setL(state_t f) { L = f; }
  inline void setPbar(state_t f) { pbar = f; }
#ifdef RESCALE
  inline state_t getScale() const { return (scale); }
#endif
  // make variables public so they can be exposed to e.g. pybind11
  // ------ variables -------
  state_t L{0};
  state_t pbar{0};
#ifdef RESCALE
  state_t scale{1.0};
#endif
private:
};

template <bool WithActivity = true>
class State : public BaseState
{
};

template <>
class State<false> : public BaseState
{
public:
  explicit State(state_t L_a = 0, state_t pbar_a = 0) : BaseState(L_a, pbar_a) {}
};

template <>
class State<true> : public BaseState
{
public:
  explicit State(state_t L_a = 0, state_t pbar_a = 0, uint8_t npa = 0, uint16_t neiq = 0)
  : BaseState(L_a, pbar_a), num_pix_active(npa), num_events_in_queue(neiq)
  {
  }

  inline uint16_t getNumEventsInQueue() const { return (num_events_in_queue); }
  inline uint8_t getNumPixActive() const { return (num_pix_active); }

  inline bool isActive() const { return (num_events_in_queue != 0); }
  inline void incNumPixActive() { num_pix_active++; }
  inline void decNumPixActive() { num_pix_active--; }
  inline void incNumEventsInQueue() { num_events_in_queue++; }
  inline void decNumEventsInQueue() { num_events_in_queue--; }
  // ------ variables -------
  uint8_t num_pix_active{0};
  uint16_t num_events_in_queue{0};
  static constexpr int max_num_active() { return (255); };  // 8 bit
};

}  // namespace fibar_lib
#endif  // FIBAR_LIB_STATE_HPP
