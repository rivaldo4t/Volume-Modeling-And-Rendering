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
	void stampWithDisplacement(lux::SField s, NoiseParams& param);
	void levelSet(Triangles& triangles);
	virtual double eval(lux::Vector p) const;
	void writeGrid(std::string fileName);
	void readGrid(std::string fileName);
};

// dirty fixes for grid scale, translate and union of large number of grids
// filling in data to be able to save chained grid operations onto file
// using f as base grid
// TODO: generalize
// eval not overloaded
class GridScale : public Grid
{
	//std::shared_ptr<Grid> f;
	//double s;
public:
	GridScale(std::shared_ptr<Grid> _f, double _s) : Grid(*_f.get())//f(_f), s(_s) {}
	{
#if 1
		gridData.resize(Nx * Ny * Nz, 0);
#pragma omp parallel for
		for (int i = 0; i < Nx; ++i)
		{
			for (int j = 0; j < Ny; ++j)
			{
				for (int k = 0; k < Nz; ++k)
				{
					lux::Vector p = llc + lux::Vector(double(i) * delta_grid, double(j) * delta_grid, double(k) * delta_grid);
					gridData[getIndex(i, j, k)] = _f->eval(p / _s);
				}
			}
		}
#endif
	}

#if 0
	double eval(lux::Vector p) const
	{
		return f->eval(p / s);
	}
#endif
};

class GridTranslate : public Grid
{
	//std::shared_ptr<Grid> f;
	//lux::Vector x;
public:
	GridTranslate(std::shared_ptr<Grid> _f, lux::Vector _x) : Grid(*_f.get())//f(_f), x(_x) {}
	{
#if 1
		gridData.resize(Nx * Ny * Nz, 0);
#pragma omp parallel for
		for (int i = 0; i < Nx; ++i)
		{
			for (int j = 0; j < Ny; ++j)
			{
				for (int k = 0; k < Nz; ++k)
				{
					lux::Vector p = llc + lux::Vector(double(i) * delta_grid, double(j) * delta_grid, double(k) * delta_grid);
					gridData[getIndex(i, j, k)] = _f->eval(p - _x);
				}
			}
		}
#endif
	}

#if 0
	double eval(lux::Vector p) const
	{
		return f->eval(p - x);
	}
#endif
};

class GridUnion : public Grid
{
	//std::shared_ptr<Grid> f, g;
public:
	GridUnion(std::shared_ptr<Grid> _f, std::shared_ptr<Grid> _g) : Grid(*_f.get())//, f(_f), g(_g) 
	{
#if 1
		gridData.resize(Nx * Ny * Nz, 0);
#pragma omp parallel for
		for (int i = 0; i < Nx; ++i)
		{
			for (int j = 0; j < Ny; ++j)
			{
				for (int k = 0; k < Nz; ++k)
				{
					lux::Vector p = llc + lux::Vector(double(i) * delta_grid, double(j) * delta_grid, double(k) * delta_grid);
					gridData[getIndex(i, j, k)] = std::max(_f->eval(p), _g->eval(p));
				}
			}
		}
		//g.reset();
#endif
	}

#if 0
	double eval(lux::Vector p) const
	{
		return std::max(f->eval(p), g->eval(p));
	}
#endif
};