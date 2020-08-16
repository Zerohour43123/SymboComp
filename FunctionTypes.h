#pragma once

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

