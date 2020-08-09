#include <iostream>
#include <math.h>
#include <string>
#include <array>
constexpr double PI = 3.14159265358979;
constexpr double NaN = ((static_cast<double>(INFINITY)) * 0);

static int funcdeletions = 0;

struct Vector2D {

	double x, y, r, theta;

	Vector2D(double x, double y) : x(x), y(y)  {
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
	std::string toString () {
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

struct FuncBase {
	virtual double evaluate(double x) const = 0;
	virtual FuncBase* differentiate() = 0;
	virtual FuncBase* antiDifferentiate() = 0;
	virtual FuncBase* Copy() = 0;
	virtual std::string toString() const = 0;
	virtual ~FuncBase() = 0;
};

FuncBase::~FuncBase() {}

struct Function : public FuncBase{
private:
	struct Monomial : public FuncBase {
		double coefficient, power;

		Monomial(double c, double p) : coefficient(c), power(p) {}

		~Monomial(){}

		Monomial(const Monomial& other) :coefficient(other.coefficient), power(other.power) {}

		FuncBase* Copy() {
			return new Monomial(*this);
		}

		std::string toString() const {
			if (power == 0 && coefficient == 0)
				return "0";
			else if (power == 0)
				return std::to_string(coefficient);
			else if (coefficient == 0)
				return "";
			else if (power == 1 && coefficient == 1)
				return "x";
			else if (power == 1)
				return std::to_string(coefficient) + "x";
			else if (coefficient == 1)
				return  "x^" + std::to_string(power);
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

		FuncBase* differentiate() {
			return new Monomial(coefficient * power, power - 1);
		}

		FuncBase* antiDifferentiate() {
			return new Monomial(coefficient / (power+1), power + 1);
		}

	};

	struct Exponential : public FuncBase {
		double coefficient, base;
		Exponential(double c, double b) : coefficient(c), base(b) {}

		~Exponential() {}

		Exponential(const Exponential& other) : coefficient(other.coefficient), base(other.base) {}

		FuncBase* Copy() {
			return new Exponential(*this);
		}

		std::string toString() const {
			if (coefficient == 1)
				return std::to_string(base) + "^x";
			else
				return std::to_string(coefficient) + "*" + std::to_string(base) + "^x";
		}

		double evaluate(double x) const {
			return coefficient * pow(base, x);
		}
		FuncBase* differentiate() {
			return new Exponential(coefficient * log(base), base);
		}

		FuncBase* antiDifferentiate() {
			return new Exponential(coefficient / log(base), base);
		}
	};

	struct Logarithm : public FuncBase {
		double coefficient;
		Function* fptr1 = 0;
		Function* fptr2 = 0;
		Function* tempptr = 0;
		bool antiDifferentiated = false;

		Logarithm(double c) : coefficient(c) {
			
		}

		~Logarithm() {
			if (antiDifferentiated) {
				fptr1 = 0;
				fptr2 = 0;
				tempptr = 0;
			}
		}
		
		Logarithm(const Logarithm& other) : coefficient(other.coefficient){
			antiDifferentiated = other.antiDifferentiated;
			if (fptr1) {
				fptr1 = new Function(*(other.fptr1));
				fptr2 = new Function(*(other.fptr2));
				tempptr = new Function(*(other.tempptr));
			}
		}

		FuncBase* Copy() {
			return new Logarithm(*this);
		}

		std::string toString() const {
			if (coefficient == 1)
				return "ln(x)";
			else
				return std::to_string(coefficient) + "*ln(x)";
		}

		double evaluate(double x) const {
			return coefficient * log(x);
		}

		FuncBase* differentiate() {
			return new Monomial(1, -1);
		}

		FuncBase* antiDifferentiate() {
			if (!antiDifferentiated) {
				antiDifferentiated = true;
				fptr1 = new Function();
				fptr1->addMode = false;
				fptr1->addMonomial(1, 1);
				fptr1->addLogarithm(1);
				fptr2 = new Function();
				fptr2->addMonomial(-1, 1);
				tempptr = new Function();
				tempptr->addFunc(fptr1);
				tempptr->addFunc(fptr2);
				return tempptr;
			}
			else
				return new Function(*tempptr);
		}
	};

	struct Trigonmetric : public FuncBase {
		double coefficient;
		char type;
		Trigonmetric(double c, char t) : coefficient(c), type(t) {
			if ((type != 's') && (type != 'c'))
				throw std::exception("Value of 'char type' must be either 's', 'c',!");
		}

		~Trigonmetric() {}

		Trigonmetric(const Trigonmetric& other) : coefficient(other.coefficient), type(other.type){}

		FuncBase* Copy() {
			return new Trigonmetric(*this);
		}

		std::string toString() const {
			if (coefficient == 1) {
				if (type == 's')
					return "sin(x)";
				else 
					return "cos(x)";
			}
			else
				if (type == 's')
					return std::to_string(coefficient) + "*sin(x)";
				else 
					return std::to_string(coefficient) + "*cos(x)";
		}



		double evaluate(double x) const {
			if (type == 's')
				return coefficient * sin(x);
			else
				return coefficient * cos(x);
		}
		
		FuncBase* differentiate() {
			if (type == 's')
				return new Trigonmetric(1, 'c');
			else
				return new Trigonmetric(-1, 's');
		}

		FuncBase* antiDifferentiate() {
			if (type == 's')
				return new Trigonmetric(-1, 'c');
			else
				return new Trigonmetric(1, 's');
		}
	};

	struct AbsoluteValue : public FuncBase {
		double coefficient;
		bool differentiated = false;
		bool antiDifferentiated = false;
		Function* fptr1 = 0;
		Function* fptr2 = 0;

		AbsoluteValue(double c) : coefficient(c) {}

		~AbsoluteValue() {
			if (differentiated)
				fptr1 = 0;
			if (antiDifferentiated)
				fptr2 = 0;
		}

		AbsoluteValue(const AbsoluteValue& other) : coefficient(other.coefficient){
			
		}

		FuncBase* Copy() {
			return new AbsoluteValue(*this);
		}

		std::string toString() const {
			if (coefficient == 1)
				return "|x|";
			else
				return std::to_string(coefficient) + "*|x|";
		}

		double evaluate(double x) const {
			return coefficient * abs(x);
		}

		FuncBase* differentiate() {
			if (!differentiated) {
				differentiated = true;
				fptr1 = new Function();
				fptr1->addMode = false;
				fptr1->addMonomial(1, -1);
				fptr1->addAbsoluteValue(1);
				return fptr1;
			}
			else
				return new Function(*fptr1);
		}

		FuncBase* antiDifferentiate() {
			if (!antiDifferentiated) {
				antiDifferentiated = true;
				fptr2 = new Function();
				fptr2->addMode = false;
				fptr2->addMonomial(.5, 1);
				fptr2->addAbsoluteValue(1);
				return fptr2;
			}
			else
				return new Function(*fptr2);
		}
	};

	//Composite Function is incomplete

	struct CompositeSubFunction : public FuncBase { 
		FuncBase* outer = 0;
		FuncBase* inner = 0;

		CompositeSubFunction(FuncBase* o, FuncBase* i) : outer(o), inner(i) {}

		~CompositeSubFunction () {
			//May or may not be useful later
			delete outer;
			delete inner;
			outer = 0;
			inner = 0;
		}

		CompositeSubFunction(const CompositeSubFunction& other) {

		}

		double evaluate(double x) const {
			return outer->evaluate(inner->evaluate(x));
		}
	};

	void addFunc(FuncBase* func) {
		if (fcount < 5) {
			functions[fcount] = func;
			fcount++;
		}
		else
			std::cout << "This function is full!\n";
	}
	
public:
	std::array<FuncBase*, 5> functions;
	unsigned int fcount = 0;
	bool addMode = true; //False is multiply


	Function() {
		for (unsigned int i = 0; i < 5; i++)
			functions[i] = 0;
	}

	~Function() {
	for (unsigned int i = 0; i < fcount; i++) 
		if (functions[i] != 0) {
			delete functions[i];
			functions[i] = 0;
		}				
	}

	Function(const Function& other)  {
		for (unsigned int i = 0; i < other.fcount; i++) {
			addFunc(other.functions[i]->Copy());
		}
		addMode = other.addMode;
	}

	FuncBase* Copy() {
		return new Function(*this);
	}

	void addMonomial(double c, double p) {
		if (fcount < 5) {
			functions[fcount] = new Monomial(c,p);
			fcount++;
		}
		else
			std::cout << "This function is full!\n";
	}

	void addExponential(double c, double b) {
		if (fcount < 5) {
			functions[fcount] = new Exponential(c, b);
			fcount++;
		}
		else
			std::cout << "This function is full!\n";
	}

	void addLogarithm(double c) {
		if (fcount < 5) {
			functions[fcount] = new Logarithm(c);
			fcount++;
		}
		else
			std::cout << "This function is full!\n";
	}

	void addTrigonmetric(double c, char t) {
		if (fcount < 5) {
			functions[fcount] = new Trigonmetric(c,t);
			fcount++;
		}
		else
			std::cout << "This function is full!\n";
	}

	void addAbsoluteValue(double c){
		if (fcount < 5) {
			functions[fcount] = new AbsoluteValue(c);
			fcount++;
		}
		else
			std::cout << "This function is full!\n";
	}

	void addComposite(){}

	double evaluate(double x) const {
		if (addMode) {
			double total = 0;
			for (unsigned int i = 0; i < fcount; i++) 
				total += functions[i]->evaluate(x);
			return total;
		}
		else {
			double total = 1;
			for (unsigned int i = 0; i < fcount; i++) 
				total *= functions[i]->evaluate(x);
			return total;
		}
		
	}

	std::string toString() const {
		std::string output;
		if (addMode) 
			for (unsigned int i = 0; i < fcount; i++) {
				if (i != 0)
					output += "+" + functions[i]->toString();
				else
					output += functions[i]->toString();
			}
		else
			for (unsigned int i = 0; i < fcount; i++) {
				if (i != 0)
					output += "*" + functions[i]->toString();
				else
					output += functions[i]->toString();
			}
		return output;
	}

	Function* differentiate() {
		Function* derivative = new Function;
		if (addMode) {
			for(unsigned int i = 0; i < fcount; i++)
				derivative->addFunc(functions[i]->differentiate());
			return derivative;
		}
		else {
				Function* dgroupsptr = new Function[fcount];
				Function* tempptr = dgroupsptr;
				for (unsigned int i = 0; i < fcount; i++) {
					for (unsigned int j = 0; j < fcount; j++) {
						if (i == j) 
							dgroupsptr->addFunc(functions[j]->differentiate());
							
						else 
							dgroupsptr->addFunc(functions[j]->Copy());

					}
					dgroupsptr->addMode = false;
					derivative->addFunc(dgroupsptr->Copy());
					dgroupsptr++;
				}
				dgroupsptr = tempptr;
				delete [] dgroupsptr;
				dgroupsptr = 0;
				return derivative;
		}
	}

	//ANTIDIFFERENTIATION NOT FULLY SUPPORTED

	Function* antiDifferentiate(){
		Function* antiDerivative = new Function();
		if (addMode) {
			for (unsigned int i = 0; i < fcount; i++)
				antiDerivative->addFunc(functions[i]->antiDifferentiate());
			return antiDerivative;
		}
		else {
			throw std::exception();
		}
	}

};

struct Calculus {

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



