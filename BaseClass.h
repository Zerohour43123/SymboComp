#pragma once

struct FuncBase {
	virtual double evaluate(double x) const = 0;
	virtual FuncBase* differentiate() const = 0;
	virtual FuncBase* antiDifferentiate() const = 0;
	virtual FuncBase* Copy() = 0;
	virtual std::string toString() const = 0;
	virtual ~FuncBase() = 0;
};

inline FuncBase::~FuncBase() {}