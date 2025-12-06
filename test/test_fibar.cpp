// -*-c++-*--------------------------------------------------------------------
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

#include <gtest/gtest.h>

#include <fibar_lib/image_reconstructor.hpp>

static const int WIDTH = 640;
static const int HEIGHT = 480;

std::shared_ptr<fibar_lib::ImageReconstructor<true, 2>> createReconstructor(
  const int w, const int h, const int T_cut, const float fill_ratio)
{
  auto reconstructor = std::make_shared<fibar_lib::ImageReconstructor<true, 2>>();
  reconstructor->initialize(w, h, T_cut, fill_ratio);
  EXPECT_EQ(w, reconstructor->getWidth());
  EXPECT_EQ(h, reconstructor->getHeight());
  EXPECT_EQ(2000, reconstructor->getEventWindowSize());
  EXPECT_EQ(0, reconstructor->getCurrentQueueSize());
  EXPECT_EQ(w * h, reconstructor->getState().size());
  EXPECT_EQ(-1.0, reconstructor->getCurrentFillRatio());  // no pix occupied yet
  return (reconstructor);
}

void test_response(const uint32_t T_cut, const float pbar, const float L)
{
  auto rec = createReconstructor(WIDTH, HEIGHT, T_cut, 0.5f);
  const uint32_t t = 0;
  const uint16_t x = WIDTH / 2;
  const uint16_t y = HEIGHT / 2;
  rec->event(t, x, y, 1);  // ON event in center of image
  const auto & state = rec->getState();
  const auto & s = state[y * WIDTH + x];
  EXPECT_NEAR(s.getPbar(), pbar, 1e-4);
  EXPECT_NEAR(s.getL(), L, 1e-4);
}

TEST(fibar_lib, test_single_updates)
{
  // for large T_cut, expect pbar close to 0, L close to 1
  test_response(10000, 0.00062812, 0.999686);
  // for small T_cut, expect pbar close to 2, L close to 1/2
  test_response(2, 2.0, 0.5857865);
}

int main(int argc, char ** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
