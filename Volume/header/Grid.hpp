#pragma once
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>
#include <memory>
#include "ScalarField.hpp"
#include "Triangle.hpp"
#include "Fspn.hpp"

class Grid
{
 protected:
	unsigned int Nx, Ny, Nz;
	double delta_grid = 0.01;
	double defaultVal = 0.0;
	lux::Vector llc, urc;
	std::vector<float> gridData;
 public:
	Grid() { Nx = 0; Ny = 0; Nz = 0; }
	Grid(lux::Vector o, unsigned int x, unsigned int y, unsigned int z, double delta, 
		std::vector<float>& data = std::vector<float>())
		: llc(o), Nx(x), Ny(y), Nz(z), delta_grid(delta)
	{ 
		gridData = std::move(data); 
		urc = { llc.X() + (Nx - 1) * delta_grid, 
				llc.Y() + (Ny - 1) * delta_grid, 
				llc.Z() + (Nz - 1) * delta_grid };
	}

	unsigned int getIndex(unsigned int i, unsigned int j, unsigned int k) const;
	bool withinGrid(lux::Vector p) const;
	void stamp(lux::SField s);
	void stampWithDisplacement(lux::SField s);
	void levelSet(Triangles& triangles);
	virtual double eval(lux::Vector p) const;
	void writeGrid(std::string fileName);
	void readGrid(std::string fileName);
};

class GridScale : public Grid
{
	std::shared_ptr<Grid> f;
	double s;
public:
	GridScale(std::shared_ptr<Grid> _f, double _s) : f(_f), s(_s) {}
	double eval(lux::Vector p) const
	{
		return f->eval(p / s);
	}
};

class GridTranslate : public Grid
{
	std::shared_ptr<Grid> f;
	lux::Vector x;
public:
	GridTranslate(std::shared_ptr<Grid> _f, lux::Vector _x) : f(_f), x(_x) {}
	double eval(lux::Vector p) const
	{
		return f->eval(p - x);
	}
};

class GridUnion : public Grid
{
	std::shared_ptr<Grid> f, g;
public:
	GridUnion(std::shared_ptr<Grid> _f, std::shared_ptr<Grid> _g) : f(_f), g(_g) {}

	double eval(lux::Vector p) const
	{
		return std::max(f->eval(p), g->eval(p));
	}
};