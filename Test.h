#include <iostream>
#include <math.h>
#include <string>
#include <array>

struct Monomial : public FuncBase {
	double coefficient, power;

	Monomial(double c, double p) : coefficient(c), power(p) {}

	~Monomial() {}

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
			return INFINITY;
		else if (power == 0)
			return coefficient;
		else
			return coefficient * pow(x, power);
	}

	FuncBase* differentiate() const {
		if (power == 0)
			return new Monomial(0, 0);
		else
			return new Monomial(coefficient * power, power - 1);
	}

	FuncBase* antiDifferentiate() const {
		if (power == -1)
			return new Composite(new Logarithm(1), new AbsoluteValue(1));
		else
			return new Monomial(coefficient / (power + 1), power + 1);
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
	FuncBase* differentiate() const {
		return new Exponential(coefficient * log(base), base);
	}

	FuncBase* antiDifferentiate() const {
		return new Exponential(coefficient / log(base), base);
	}
};

struct Logarithm : public FuncBase {
	double coefficient;

	Logarithm(double c) : coefficient(c) {

	}

	~Logarithm() {}

	Logarithm(const Logarithm& other) : coefficient(other.coefficient) {}

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

	FuncBase* differentiate() const {
		return new Monomial(1, -1);
	}

	FuncBase* antiDifferentiate() const {
		Function* antiDerivative = new Function();
		Function part1;
		part1.numeratorAddMode = false;
		part1.addNumerMonomial(1, 1);
		part1.addNumerLogarithm(1);
		Function part2;
		part2.addNumerMonomial(-1, 1);
		antiDerivative->addNumerFunc(part1.Copy());
		antiDerivative->addNumerFunc(part2.Copy());
		return antiDerivative;

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

	Trigonmetric(const Trigonmetric& other) : coefficient(other.coefficient), type(other.type) {}

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

	FuncBase* differentiate() const {
		if (type == 's')
			return new Trigonmetric(1, 'c');
		else
			return new Trigonmetric(-1, 's');
	}

	FuncBase* antiDifferentiate() const {
		if (type == 's')
			return new Trigonmetric(-1, 'c');
		else
			return new Trigonmetric(1, 's');
	}
};

struct AbsoluteValue : public FuncBase {
	double coefficient;

	AbsoluteValue(double c) : coefficient(c) {}

	~AbsoluteValue() {}

	AbsoluteValue(const AbsoluteValue& other) : coefficient(other.coefficient) {}

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

	FuncBase* differentiate() const {
		Function* derivative = new Function();
		derivative->numeratorAddMode = false;
		derivative->addNumerMonomial(1, -1);
		derivative->addNumerAbsoluteValue(1);
		return derivative;
	}

	FuncBase* antiDifferentiate() const {
		Function* antiDerivative = new Function();
		antiDerivative->numeratorAddMode = false;
		antiDerivative->addNumerMonomial(.5, 1);
		antiDerivative->addNumerAbsoluteValue(1);
		return antiDerivative;
	}
};

struct Composite : public FuncBase {
	FuncBase* outer = 0;
	FuncBase* inner = 0;

	Composite(FuncBase* o, FuncBase* i) : outer(o), inner(i) {
	}

	~Composite() {
		delete outer;
		outer = 0;
		delete inner;
		inner = 0;
	}

	std::string toString() const {
		int position = outer->toString().find("x");
		return outer->toString().substr(0, position) + "(" + inner->toString() + ")" + outer->toString().substr(position + 1);
	}

	Composite(const Composite& other) {
		outer = other.outer->Copy();
		inner = other.outer->Copy();
	}

	FuncBase* Copy() {
		return new Composite(*this);
	}

	double evaluate(double x) const {
		return outer->evaluate(inner->evaluate(x));
	}

	FuncBase* differentiate() const {
		Function* derivative = new Function();
		derivative->numeratorAddMode = false;
		derivative->addNumerComposite(outer->differentiate(), inner->Copy());
		derivative->addNumerFunc(inner->differentiate());
		return derivative;
	}

	FuncBase* antiDifferentiate() const {
		FuncBase* tempptr = inner->differentiate();
		if (tempptr->evaluate(49.3418) == tempptr->evaluate(374.92687)) {
			Function* antiDerivative = new Function();
			antiDerivative->numeratorAddMode = false;
			antiDerivative->addNumerComposite(outer->antiDifferentiate(), inner->Copy());
			antiDerivative->addNumerComposite(new Monomial(1, -1), inner->differentiate());
			delete tempptr;
			return antiDerivative;
		}
		else
			throw std::exception();
	}
};

#pragma once
#include <iostream>
#include <math.h>
#include <string>
#include <array>
#include "Calculus.h"
#include "Vector.h"
#include "FunctionTypes.h"

constexpr double PI = 3.1415926535897932;
constexpr double e = 2.7182818284590451;
constexpr double NaN = ((static_cast<double>(INFINITY)) * 0);



struct FuncBase {
	virtual double evaluate(double x) const = 0;
	virtual FuncBase* differentiate() const = 0;
	virtual FuncBase* antiDifferentiate() const = 0;
	virtual FuncBase* Copy() = 0;
	virtual std::string toString() const = 0;
	virtual ~FuncBase() = 0;
};

inline FuncBase::~FuncBase() {}



struct Function : public FuncBase {



	std::array<FuncBase*, 5> numeratorFunctions;
	std::array<FuncBase*, 5> denominatorFunctions;
	unsigned int numfcount = 0;
	unsigned int denomfcount = 0;
	bool numeratorAddMode = true; //False is multiply
	bool denominatorAddMode = true;     //False is multiply

	Function() {
		for (unsigned int i = 0; i < 5; i++) {
			numeratorFunctions[i] = 0;
			denominatorFunctions[i] = 0;
		}
	}

	~Function() {
		for (unsigned int i = 0; i < numfcount; i++) {
			if (numeratorFunctions[i] != 0) {
				delete numeratorFunctions[i];
				numeratorFunctions[i] = 0;
			}
		}
		for (unsigned int i = 0; i < denomfcount; i++) {
			if (denominatorFunctions[i] != 0) {
				delete denominatorFunctions[i];
				denominatorFunctions[i] = 0;
			}
		}
	}

	Function(const Function& other) {
		for (unsigned int i = 0; i < other.numfcount; i++) {
			addNumerFunc(other.numeratorFunctions[i]->Copy());
		}
		numeratorAddMode = other.numeratorAddMode;
		for (unsigned int i = 0; i < other.denomfcount; i++) {
			addDenomFunc(other.denominatorFunctions[i]->Copy());
		}
		denominatorAddMode = other.denominatorAddMode;
	}

	FuncBase* Copy() {
		return new Function(*this);
	}

	void addNumerMonomial(double c, double p) {
		if (numfcount < 5) {
			numeratorFunctions[numfcount] = new Monomial(c, p);
			numfcount++;
		}
		else
			std::cout << "This function is full!\n";
	}

	void addNumerExponential(double c, double b) {
		if (numfcount < 5) {
			numeratorFunctions[numfcount] = new Exponential(c, b);
			numfcount++;
		}
		else
			std::cout << "This function is full!\n";
	}

	void addNumerLogarithm(double c) {
		if (numfcount < 5) {
			numeratorFunctions[numfcount] = new Logarithm(c);
			numfcount++;
		}
		else
			std::cout << "This function is full!\n";
	}

	void addNumerTrigonmetric(double c, char t) {
		if (numfcount < 5) {
			numeratorFunctions[numfcount] = new Trigonmetric(c, t);
			numfcount++;
		}
		else
			std::cout << "This function is full!\n";
	}

	void addNumerAbsoluteValue(double c) {
		if (numfcount < 5) {
			numeratorFunctions[numfcount] = new AbsoluteValue(c);
			numfcount++;
		}
		else
			std::cout << "This function is full!\n";
	}

	void addNumerComposite(FuncBase* o, FuncBase* i) {
		if (numfcount < 5) {
			numeratorFunctions[numfcount] = new Composite(o, i); //POSSIBLY PROBLEMATIC WITH "NEW"???
			numfcount++;
		}
		else
			std::cout << "This function is full!\n";
	}

	void addNumerFunc(FuncBase* func) {
		if (numfcount < 5) {
			numeratorFunctions[numfcount] = func;
			numfcount++;
		}
		else
			std::cout << "This function is full!\n";
	}

	void addDenomMonomial(double c, double p) {
		if (denomfcount < 5) {
			denominatorFunctions[denomfcount] = new Monomial(c, p);
			denomfcount++;
		}
		else
			std::cout << "This function is full!\n";
	}

	void addDenomExponential(double c, double b) {
		if (denomfcount < 5) {
			denominatorFunctions[denomfcount] = new Exponential(c, b);
			denomfcount++;
		}
		else
			std::cout << "This function is full!\n";
	}

	void addDenomLogarithm(double c) {
		if (denomfcount < 5) {
			denominatorFunctions[denomfcount] = new Logarithm(c);
			denomfcount++;
		}
		else
			std::cout << "This function is full!\n";
	}

	void addDenomTrigonmetric(double c, char t) {
		if (denomfcount < 5) {
			denominatorFunctions[denomfcount] = new Trigonmetric(c, t);
			denomfcount++;
		}
		else
			std::cout << "This function is full!\n";
	}

	void addDenomAbsoluteValue(double c) {
		if (denomfcount < 5) {
			denominatorFunctions[denomfcount] = new AbsoluteValue(c);
			denomfcount++;
		}
		else
			std::cout << "This function is full!\n";
	}

	void addDenomComposite(FuncBase* o, FuncBase* i) {
		if (denomfcount < 5) {
			denominatorFunctions[denomfcount] = new Composite(o, i); //POSSIBLY PROBLEMATIC WITH "NEW"???
			denomfcount++;
		}
		else
			std::cout << "This function is full!\n";
	}

	void addDenomFunc(FuncBase* func) {
		if (denomfcount < 5) {
			denominatorFunctions[denomfcount] = func;
			denomfcount++;
		}
		else
			std::cout << "This function is full!\n";
	}

	void movetoNumerator(int position) {
		if (typeid(*(denominatorFunctions[position])).name()[17] == 'M') {
			((Monomial*)denominatorFunctions[position])->power *= -1;
			addNumerFunc(denominatorFunctions[position]);
		}
		else
			addNumerComposite(new Monomial(1, -1), numeratorFunctions[position]);
		if (position + 1 == denomfcount) {
			denominatorFunctions[position] = 0;
			denomfcount--;
		}
		else {
			denomfcount--;
			for (unsigned int i = 0; i < denomfcount; i++)
				denominatorFunctions[i] = denominatorFunctions[i + 1];
			denominatorFunctions[position + 1] = 0;
		}
	}

	void movetoDenominator(int position) {
		if (typeid(*(numeratorFunctions[position])).name()[17] == 'M') {
			((Monomial*)numeratorFunctions[position])->power *= -1;
			addDenomFunc(numeratorFunctions[position]);
		}
		else
			addDenomComposite(new Monomial(1, -1), numeratorFunctions[position]);
		if (position == numfcount) {
			numeratorFunctions[position] = 0;
			numfcount--;
		}
		else {
			numfcount--;
			for (unsigned int i = 0; i < numfcount; i++)
				numeratorFunctions[i] = numeratorFunctions[i + 1];
			numeratorFunctions[position + 1] = 0;
		}
	}

	void simplify() {
		if (numfcount < 2 || numeratorAddMode)
			throw std::exception("\"simplify()\" requires numerator functions in multiply mode!");
		unsigned int monomialCount = 0;
		unsigned int firstMonomial = 0;
		double newPower = 0;
		double newCoefficient = 1;
		for (unsigned int i = 0; i < numfcount; i++) {
			if (typeid(*(numeratorFunctions[i])).name()[17] == 'M') {
				if (monomialCount == 0) {
					firstMonomial = i;
					monomialCount++;
				}
				else {
					newPower += (static_cast<Monomial*>(numeratorFunctions[i]))->power;
					newCoefficient *= (static_cast<Monomial*>(numeratorFunctions[i]))->coefficient;
					delete numeratorFunctions[i];
					if (i != (numfcount - 1))
						numeratorFunctions[i] = numeratorFunctions[i + 1];
					else
						numeratorFunctions[i] = 0;
					monomialCount++;
				}
			}
		}
		static_cast<Monomial*>(numeratorFunctions[firstMonomial])->power += newPower;
		static_cast<Monomial*>(numeratorFunctions[firstMonomial])->coefficient *= newCoefficient;
		numfcount -= (monomialCount - 1);
	}

	double evaluate(double x) const {
		double numerTotal = 0;
		if (numeratorAddMode) {
			for (unsigned int i = 0; i < numfcount; i++)
				numerTotal += numeratorFunctions[i]->evaluate(x);
		}
		else {
			numerTotal = 1;
			for (unsigned int i = 0; i < numfcount; i++)
				numerTotal *= numeratorFunctions[i]->evaluate(x);
		}
		if (denomfcount != 0) {
			double denomTotal = 0;
			if (denominatorAddMode) {
				for (unsigned int i = 0; i < denomfcount; i++)
					denomTotal += denominatorFunctions[i]->evaluate(x);
			}
			else {
				double numerTotal = 1;
				for (unsigned int i = 0; i < denomfcount; i++)
					denomTotal *= denominatorFunctions[i]->evaluate(x);
			}
			return numerTotal / denomTotal;
		}
		else
			return numerTotal;
	}

	std::string toString() const {
		std::string output;
		if (numeratorAddMode)
			for (unsigned int i = 0; i < numfcount; i++) {
				if (i != 0)
					output += "+" + numeratorFunctions[i]->toString();
				else
					output += numeratorFunctions[i]->toString();
			}
		else
			for (unsigned int i = 0; i < numfcount; i++) {
				if (i != 0)
					output += "*" + numeratorFunctions[i]->toString();
				else
					output += numeratorFunctions[i]->toString();
			}
		return output;
	}

	Function* differentiate() const {
		Function* derivative = new Function;
		if (denomfcount != 0)
			throw std::exception("Cannot have functions in denominator whilst differentiating! s(for now)");
		if (numeratorAddMode) {
			for (unsigned int i = 0; i < numfcount; i++)
				derivative->addNumerFunc(numeratorFunctions[i]->differentiate());
			return derivative;
		}
		else {
			std::array<Function*, 5> dgroups;
			for (unsigned int i = 0; i < numfcount; i++) {
				(dgroups)[i] = new Function();
			}
			for (unsigned int i = 0; i < numfcount; i++) {
				for (unsigned int j = 0; j < numfcount; j++) {
					if (i == j)
						(dgroups)[i]->addNumerFunc(numeratorFunctions[j]->differentiate());
					else
						(dgroups)[i]->addNumerFunc(numeratorFunctions[j]->Copy());
				}
				(dgroups)[i]->numeratorAddMode = false;
				derivative->addNumerFunc((dgroups)[i]);
			}
			return derivative;
		}
	}

	//ANTIDIFFERENTIATION NOT FULLY SUPPORTED

	Function* antiDifferentiate() const {
		if (numfcount == 1 && typeid(*(numeratorFunctions[0])).name()[17] == 'L')
			return static_cast<Function*>(numeratorFunctions[0]->antiDifferentiate());
		if (numeratorAddMode) {
			Function* antiDerivative = new Function();
			for (unsigned int i = 0; i < numfcount; i++)
				antiDerivative->addNumerFunc(numeratorFunctions[i]->antiDifferentiate());
			return antiDerivative;
		}
		else {
			throw std::exception();
		}
	}

};






