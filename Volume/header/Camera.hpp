#pragma once
#include "Vector.hpp"

class Camera
{
  public:
    Camera();
   ~Camera(){}

    void setEyeViewUp( const lux::Vector& eye, const lux::Vector& view, const lux::Vector& up );
    const lux::Vector& eye() const  { return position; }
    const lux::Vector& view() const { return axis_view; }
    const lux::Vector& up() const   { return axis_up; }
	const lux::Vector& right() const { return axis_right; }

    const lux::Vector view( const double x, const double y );

    void setFov( const double fov );
    const double& fov() const { return FOV; }
	const double& htanFOV() const { return htanfov; }
	const double& vtanFOV() const { return vtanfov; }

    void setAspectRatio( const double ar );
    const double& aspectRatio() const { return aspect_ratio; }

    void setNearPlane( const double n ){ near = n; }
    const double& nearPlane() const { return near; }

    void setFarPlane( const double n ){ far = n; }
    const double& farPlane() const { return far; }

  private:
    double FOV, aspect_ratio;
    double htanfov, vtanfov;
    double near, far;
    lux::Vector position;
    lux::Vector axis_right, axis_up, axis_view;
};