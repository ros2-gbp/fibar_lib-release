# FIBAR: filter based image reconstruction from events

This repository holds plain C++ code for event image reconstruction by
means of a temporal and spatial filtering algorithm
[described here](https://arxiv.org/abs/2510.20071). It has no dependencies
on ROS.

For usage examples, see the corresponding
[ROS image reconstruction repo](https://github.com/ros-event-camera/event_image_reconstruction_fibar).

## Supported platforms

Currently tested for long-term Ubuntu 22.04 and later.

## How to build

Set the following shell variables:

```bash
repo=fibar_lib
url=https://github.com/ros-event-camera/${repo}.git
```

and follow the generic ROS package build
[instructions here](https://github.com/ros-misc-utilities/.github/blob/master/docs/build_ros_repository.md).
Since this package is a ROS-free cmake package, you can also just use the standard cmake build commands.

## How to run tests

Before committing changes to the repo, run the tests as follows at the top of the ROS workspace.

```bash
rm -rf build_test; mkdir build_test; cd build_test; cmake ../src/fibar_lib/ -DBUILD_TESTING=ON; make; make test; make format_check; cd ..
```

## License

This software is issued under the Apache License Version 2.0.
