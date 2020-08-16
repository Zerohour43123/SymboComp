#pragma once

struct Vector2D {

	double x, y, r, theta;

	Vector2D(double x, double y) : x(x), y(y) {
		r = sqrt(pow(x, 2) + pow(y, 2));
		theta = (180 / PI) * atan2(y, x);
	}
	Vector2D operator+(const Vector2D& input) const {
		return Vector2D(x + input.x, y + input.y);
	}
	Vector2D operator-(const Vector2D& input) const {
		return Vector2D(x - input.x, y - input.y);
	}
	Vector2D operator*(double scalar) const {
		return Vector2D(x * scalar, y * scalar);
	}
	std::string toString() {
		return "Component x: " + std::to_string(x) + ", Component y: " + std::to_string(y) + "\nMagnitude: " + std::to_string(r) + ", Angle theta: " + std::to_string(theta) + "\n";
	}
	double dotProduct(const Vector2D& input) const {
		return (x * input.x) + (y * input.y);
	}
	double angleBetween(const Vector2D& input) const {
		return (180 / PI) * acos(dotProduct(input) / (r * input.r));
	}
	Vector2D vectorProject(const Vector2D& input) const {
		return Vector2D(((*this) * (dotProduct(input) / pow(r, 2))).x, ((*this) * (dotProduct(input) / pow(r, 2))).y);
	}
};

struct Vector3D {

	double x, y, z, r, theta, phi;

	Vector3D(double x, double y, double z) : x(x), y(y), z(z) {
		r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
		phi = (180 / PI) * atan2(y, x);
		theta = (180 / PI) * acos(z / r);
	}

	Vector3D(const Vector3D& other) {
		x = other.x;
		y = other.y;
		z = other.z;
		r = other.r;
		theta = other.theta;
		phi = other.phi;
		printf("Copied!\n");
	}

	std::string toString() {
		return "Component x: " + std::to_string(x) + ", Component y: " + std::to_string(y) + ", Component z: " + std::to_string(z) + "\nMagnitude: " + std::to_string(r) + ", Angle theta: " + std::to_string(theta) + ", Angle phi: " + std::to_string(phi) + "\n";
	}
	Vector3D operator+(const Vector3D& input) const {
		return Vector3D(x + input.x, y + input.y, z + input.z);
	}
	Vector3D operator-(const Vector3D& input) const {
		return Vector3D(x - input.x, y - input.y, z - input.z);
	}
	Vector3D operator*(double scalar) const {
		return Vector3D(x * scalar, y * scalar, z * scalar);
	}
	double dotProduct(const Vector3D& input) const {
		return (x * input.x) + (y * input.y) + (z * input.z);
	}
	double angleBetween(const Vector3D& input) const {
		return (180 / PI) * acos(dotProduct(input) / (r * input.r));
	}
	Vector3D vectorProject(const Vector3D& input) const {
		return Vector3D(((*this) * (dotProduct(input) / pow(r, 2))).x, ((*this) * (dotProduct(input) / pow(r, 2))).y, ((*this) * (dotProduct(input) / pow(r, 2))).z);
	}
	Vector3D crossProduct(const Vector3D& input) const {
		return Vector3D(y * input.z - z * input.y, z * input.x - x * input.z, x * input.y - y * input.x);
	}
};