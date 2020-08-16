#pragma once

struct Function : public FuncBase {

#include "FunctionTypes.h"

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

