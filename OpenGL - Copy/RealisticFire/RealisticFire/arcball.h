#ifndef _ARCBALL_H
#define _ARCBALL_H

#include "core.h"

typedef float GLMatrix[16];

class Arcball
{
private:

	//window size
	int _winX, _winY;

	//for zooming
	float _zoomRate;
	float _startZoomX, _startZoomY;

	//for dragging (points in 3d scenes)
	float _dragRate;
	float _startDragX, _startDragY;

	Eigen::Vector3f _viewDir;
	Eigen::Vector3f _upDir;
	Eigen::Vector3f _rightDir;

	GLMatrix _startMatrix;
	Eigen::Vector3d _startRotationVector;
	Eigen::Vector3d _currentRotationVector;


	bool _isZooming;
	bool _isRotating;
	bool _isDragging;
	float _ballRadius;

	//tool functions
	Eigen::Vector3d ConvertXY(int x, int y);
	void ApplyRotationMatrix();
public:
	Arcball();

	void SetWidthHeight(int w, int h);
	void SetRadius(float newRadius);

	//rotation
	void StartRotation(int x, int y);
	void UpdateRotation(int x, int y);
	void StopRotation();

	//zooming
	void StartZooming(int x, int y);
	void UpdateZooming(int x, int y);
	void StopZooming();

	//dragging
	void StartDragging(int x, int y);
	Eigen::Vector3f UpdateDragging(int x, int y);
	void StopDragging();


	void Reset();
};


#endif
