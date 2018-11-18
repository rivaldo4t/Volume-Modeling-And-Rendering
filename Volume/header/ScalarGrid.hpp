#pragma once
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>
#include <memory>
#include "Field.hpp"
#include "ScalarField.hpp"
#include "Triangle.hpp"
#include "Fspn.hpp"

namespace lux
{
	class ScalarGrid : public Field<float>
	{
	protected:
		unsigned int Nx, Ny, Nz;
		double delta_grid = 0.01;
		double defaultVal = 0.0;
		lux::Vector llc, urc;
		std::vector<float> gridData;
	public:
		ScalarGrid() { Nx = 0; Ny = 0; Nz = 0; }
		ScalarGrid(lux::Vector o, unsigned int x, unsigned int y, unsigned int z, double delta,
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
		void writeGrid(std::string fileName);
		void readGrid(std::string fileName);

		virtual const float eval(const Vector& p) const override;
		void stamp(lux::SField s);
		void stampWithDisplacement(lux::SField s, NoiseParams& param);
		void pyroDisplace(NoiseParams& param);
		void levelSet(Triangles& triangles);
	};

	// dirty fixes for grid scale, translate and union of large number of grids
	// filling in data to be able to save chained grid operations onto file
	// using f as base grid
	// TODO: generalize
	// eval not overloaded
	class GridScale : public ScalarGrid
	{
		//std::shared_ptr<Grid> f;
		//double s;
	public:
		GridScale(std::shared_ptr<ScalarGrid> _f, double _s) : ScalarGrid(*_f.get())//f(_f), s(_s) {}
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

	class GridTranslate : public ScalarGrid
	{
		//std::shared_ptr<Grid> f;
		//lux::Vector x;
	public:
		GridTranslate(std::shared_ptr<ScalarGrid> _f, lux::Vector _x) : ScalarGrid(*_f.get())//f(_f), x(_x) {}
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

	class GridUnion : public ScalarGrid
	{
		//std::shared_ptr<Grid> f, g;
	public:
		GridUnion(std::shared_ptr<ScalarGrid> _f, std::shared_ptr<ScalarGrid> _g) : ScalarGrid(*_f.get())//, f(_f), g(_g) 
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

	class GridRotate : public ScalarGrid
	{
		//std::shared_ptr<Grid> f, g;
	public:
		GridRotate(std::shared_ptr<ScalarGrid> _f, double angle, lux::Vector axis) : ScalarGrid(*_f.get())//, f(_f), g(_g) 
		{
#if 1
			gridData.resize(Nx * Ny * Nz, 0);
			angle = -angle * M_PI / 180.0;
#pragma omp parallel for
			for (int i = 0; i < Nx; ++i)
			{
				for (int j = 0; j < Ny; ++j)
				{
					for (int k = 0; k < Nz; ++k)
					{
						lux::Vector p = llc + lux::Vector(double(i) * delta_grid, double(j) * delta_grid, double(k) * delta_grid);
						double A = cos(angle);
						double BB = (axis * p) * (1 - A);
						double C = sin(angle);
						lux::Vector rotated_p = (A * p) + (BB * axis) + (C * (p ^ axis));
						gridData[getIndex(i, j, k)] = _f->eval(rotated_p);
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
}