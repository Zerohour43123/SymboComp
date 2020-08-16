#pragma once

struct Calculus {

	static double Integrate(const FuncBase& func, double a, double b, bool PV) {
		double sum = -4444444;
		double c = 0;
		bool improper1 = false;
		bool improper2 = false;
		bool improper3 = false;
		FuncBase* antiDerivative = func.antiDifferentiate();
		if (isinf(a) || isinf(b))
			improper1 = true;
		if (isinf(func.evaluate(a)) || isinf(func.evaluate(b)))
			improper2 = true;
		if (!improper1) {
			for (double i = a + 1; i < b; i++)
				if (isinf(func.evaluate(i))) {
					if (improper2)
						throw std::exception("Integral must not have discontinuities at the bounds and in the middle!");
					else {
						improper3 = true;
						c = i;
						break;
					}
				}
		}

		/*if (improper1 && improper3) {
			if (isinf(a) && isinf(b))
				sum = (LeftLimit(*antiDerivative, c) - Limit(*antiDerivative, a)) + (Limit(*antiDerivative,b) - RightLimit(*antiDerivative, c));
			else if(!isinf(a))
				sum = (LeftLimit(*antiDerivative, c) - antiDerivative->evaluate(a)) + (Limit(*antiDerivative, b) - RightLimit(*antiDerivative, c));
			else
				sum = (LeftLimit(*antiDerivative, c) - Limit(*antiDerivative, a)) + (antiDerivative->evaluate(b) - RightLimit(*antiDerivative, c));
		}
		else */ if (improper1) {
			if (isinf(a) && isinf(b) && PV)
				sum = (Limit(*antiDerivative, b) - Limit(*antiDerivative, a));
			else if (isinf(a) && isinf(b))
				sum = (antiDerivative->evaluate(0) - Limit(*antiDerivative, a)) + (Limit(*antiDerivative, b) - antiDerivative->evaluate(0));
			else if (!isinf(a))
				sum = (Limit(*antiDerivative, b) - antiDerivative->evaluate(a));
			else
				sum = (antiDerivative->evaluate(b) - Limit(*antiDerivative, a));
		}
		else if (improper2) {
			if (isinf(func.evaluate(a)) && isinf(func.evaluate(b)) && PV)
				sum = -3333333;
			else if (isinf(func.evaluate(a)) && isinf(func.evaluate(b))) {
				if (typeid(*(antiDerivative)).name()[7] == 'F' && static_cast<Function*>(antiDerivative)->numfcount > 1) {
					double lLimitSum = 0;
					double rLimitSum = 0;
					for (unsigned int i = 0; i < static_cast<Function*>(antiDerivative)->numfcount; i++) {
						lLimitSum += LeftLimit(*(static_cast<Function*>(antiDerivative)->numeratorFunctions[i]), b);
						rLimitSum += RightLimit(*(static_cast<Function*>(antiDerivative)->numeratorFunctions[i]), a);
					}
					sum = (antiDerivative->evaluate((b - a) / 2) - rLimitSum) + (lLimitSum - antiDerivative->evaluate((b - a) / 2));;
				}
				else
					sum = (antiDerivative->evaluate((b - a) / 2) - RightLimit(*antiDerivative, a)) + (LeftLimit(*antiDerivative, b) - antiDerivative->evaluate((b - a) / 2));
			}
			else if (!isinf(func.evaluate(a))) {
				if (typeid(*(antiDerivative)).name()[7] == 'F' && static_cast<Function*>(antiDerivative)->numfcount > 1) {
					double lLimitSum = 0;
					for (unsigned int i = 0; i < static_cast<Function*>(antiDerivative)->numfcount; i++) {
						lLimitSum += LeftLimit(*(static_cast<Function*>(antiDerivative)->numeratorFunctions[i]), b);
					}
					sum = lLimitSum - antiDerivative->evaluate(b);
				}
				else
					sum = (LeftLimit(*antiDerivative, b) - antiDerivative->evaluate(a));
			}
			else {
				if (typeid(*(antiDerivative)).name()[7] == 'F' && static_cast<Function*>(antiDerivative)->numfcount > 1) {
					double rLimitSum = 0;
					for (unsigned int i = 0; i < static_cast<Function*>(antiDerivative)->numfcount; i++) {
						rLimitSum += RightLimit(*(static_cast<Function*>(antiDerivative)->numeratorFunctions[i]), a);
					}
					sum = antiDerivative->evaluate(b) - rLimitSum;
				}
				else
					sum = (antiDerivative->evaluate(b) - RightLimit(*antiDerivative, a));
			}

		}
		else if (improper3) {
			sum = (LeftLimit(*antiDerivative, c) - antiDerivative->evaluate(a)) + (antiDerivative->evaluate(b) - RightLimit(*antiDerivative, c));
		}
		else
			sum = (antiDerivative->evaluate(b) - antiDerivative->evaluate(a));

		delete antiDerivative;
		antiDerivative = 0;
		return sum;
	}

	static double Limit(FuncBase& func, double x) {
		double evaluation = func.evaluate(x);
		if (isnan(evaluation)) {
			double position = NaN;
			Function* flhop = (Function*)(func.Copy());
			if (flhop->denomfcount == 0)
				flhop->movetoDenominator(0);
			Function toSimplify;
			toSimplify.addNumerFunc(flhop->numeratorFunctions[0]->differentiate());
			toSimplify.addDenomFunc(flhop->denominatorFunctions[0]->differentiate());
			toSimplify.numeratorAddMode = false;
			toSimplify.movetoNumerator(0);
			toSimplify.simplify();
			position = Limit(toSimplify, x);
			delete flhop;
			flhop = 0;
			return position;
		}
		else if (isinf(x) && isinf(evaluation))
			return evaluation;
		else if (!isinf(x) && isinf(evaluation)) {
			double leftval = LeftLimit(func, x);
			double rightval = RightLimit(func, x);
			if (leftval == rightval)
				return leftval;
			else
				return NaN;
		}
		else
			return evaluation;
	};

	static double LeftLimit(FuncBase& func, double x) {
		double evaluation = func.evaluate(x);
		if (isinf(x))
			throw std::runtime_error("Left limit of " + std::to_string(x) + " is not possible!");
		else if (isnan(evaluation)) {
			double position = NaN;
			Function* flhop = (Function*)(func.Copy());
			if (flhop->denomfcount == 0)
				flhop->movetoDenominator(0);
			Function toSimplify;
			toSimplify.addNumerFunc(flhop->numeratorFunctions[0]->differentiate());
			toSimplify.addDenomFunc(flhop->denominatorFunctions[0]->differentiate());
			toSimplify.numeratorAddMode = false;
			toSimplify.movetoNumerator(0);
			toSimplify.simplify();
			position = LeftLimit(toSimplify, x);
			delete flhop;
			flhop = 0;
			return position;
		}
		else if (isinf(evaluation)) {
			if (func.evaluate(x - .001) > 0)
				return INFINITY;
			else
				return -INFINITY;
		}
		else
			return evaluation;
	}

	static double RightLimit(FuncBase& func, double x) {
		double evaluation = func.evaluate(x);
		if (isinf(x))
			throw std::runtime_error("Right limit of " + std::to_string(x) + " is not possible!");
		else if (isnan(evaluation)) {
			double position = NaN;
			Function* flhop = (Function*)(func.Copy());
			if (flhop->denomfcount == 0)
				flhop->movetoDenominator(0);
			Function toSimplify;
			toSimplify.addNumerFunc(flhop->numeratorFunctions[0]->differentiate());
			toSimplify.addDenomFunc(flhop->denominatorFunctions[0]->differentiate());
			toSimplify.numeratorAddMode = false;
			toSimplify.movetoNumerator(0);
			toSimplify.simplify();
			position = RightLimit(toSimplify, x);
			delete flhop;
			flhop = 0;
			return position;
		}
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