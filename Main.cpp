#include <iostream>
#include <math.h>
#include <string>
constexpr double PI = 3.14159265358979;
constexpr double NaN = ((static_cast<double>(INFINITY)) * 0);

struct Vector2D {

	double x, y, r, theta;

	Vector2D(double x, double y) {
		this->x = x;
		this->y = y;
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
	void printVector() {
		std::cout << ("Component x: " + std::to_string(x) + ", Component y: " + std::to_string(y) + "\nMagnitude: " + std::to_string(r) + ", Angle theta: " + std::to_string(theta) + "\n");
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

	Vector3D(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
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
	void printVector() {
		std::cout << ("Component x: " + std::to_string(x) + ", Component y: " + std::to_string(y) + ", Component z: " + std::to_string(z) + "\nMagnitude: " + std::to_string(r) + ", Angle theta: " + std::to_string(theta) + ", Angle phi: " + std::to_string(phi) + "\n");
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

struct Function {
	virtual double evaluate(double x) const = 0;
	virtual Function* differentiate() const = 0;
	virtual Function* antidifferentiate() const = 0;
	virtual ~Function() = 0;
};

Function::~Function() {}

struct Monomial : public Function {
	double coefficient, power;
	Monomial(double c, double p) {
		coefficient = c;
		power = p;
	}
	~Monomial() {}

	std::string toString() {
		if (power == 0)
			return std::to_string(coefficient);
		else if (coefficient == 0)
			return "";
		else if (power == 1 && coefficient == 1)
			return "x";
		else if (power == 1)
			return std::to_string(coefficient) + "x";
		else
			return std::to_string(coefficient) + "*x^" + std::to_string(power);
	}

	double evaluate(double x) const {
		if (power < 0 && x == 0)
			return -INFINITY;
		else if (fmod(abs(power), 2) == 0 && x < 0)
			return NaN;
		else
			return coefficient * pow(x, power);
	}

	Function* differentiate() const {
		return new Monomial(coefficient * power, power - 1);
	}
	Function* antidifferentiate() const {
		return new Monomial(coefficient / power, power + 1);
	}

};

class PolyFunction : public Function {

};

class Calculus {
public:
	static double Limit(const Function& func, double x) {
		double evaluation = func.evaluate(x);
		if (isnan(evaluation))
			return nan("Complex!");
		else if (isinf(x) && isinf(evaluation))
			return evaluation;
		else if (!isinf(x) && isinf(evaluation)) {
			return 0;
		}
		else
			return evaluation;
	};

	static double LeftLimit(const Function& func, double x) {
		double evaluation = func.evaluate(x);
		if (isinf(x))
			throw std::runtime_error("Left limit of " + std::to_string(x) + " is not possible!");
		else if (isnan(evaluation))
			return nan("Complex!");
		else if (isinf(evaluation)) {
			if (func.evaluate(x - .001) > 0)
				return INFINITY;
			else
				return -INFINITY;
		}
		else
			return evaluation;
	}

	static double RightLimit(const Function& func, double x) {
		double evaluation = func.evaluate(x);
		if (isinf(x))
			throw std::runtime_error("Right limit of " + std::to_string(x) + " is not possible!");
		else if (isnan(evaluation))
			return nan("Complex!");
		else if (isinf(evaluation)) {
			if (func.evaluate(x + .001) > 0)
				return INFINITY;
			else
				return -INFINITY;
		}
		else
			return evaluation;
	}



};

int main()
{
	return 0;
}


