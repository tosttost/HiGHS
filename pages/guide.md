---
permalink: /guide/
---

## Guide

This guide gives a brief overview of what HiGHS can do when called
from an application rather than run as an
[executable](#run-executable). Detailed formal documentation setting
out the data structures and methods in the C++ interface, the
corresponding methods in the C interface, and how HiGHS can be
controlled, is [available](../docs/HighsDocumentation.pdf "HiGHS
Documentation"). The methods in the FORTRAN 90 and C# interfaces are
functionally identical to the methods in the C interface.

### Model definition

Models are either passed to HiGHS directly, read from data files, or
built within HiGHS by adding columns and rows to an empty model. When
the model is passed or built, the matrix information is communicated
in [packed vector form](https://en.wikipedia.org/wiki/Sparse_matrix),
either row-wise or column-wise. At any point in time HiGHS has a
single *incumbent model*.

#### Passing a model to HiGHS

When called from C++, a model is passed as an instance of the
`HighsModel` class (in general), or as an instance of the `HighsLp`
class if the model is a (mixed integer) linear programming problem. In
the C interface, the model is defined by passing the corresponding
scalar and array data as parameters.

#### Reading from a data file

HiGHS can read a model from files in the MPS or CPLEX LP formats. If
an MPS file contains a section corresponding to an extension not
supported by HiGHS, an error will be returned. For models in the CPLEX
LP format, HiGHS will read no further than the standard maximum line
length of 255 characters. HiGHS will return with an error if longer
lines are encountered.

#### Building a model

In a newly-created instance of HiGHS, the incumbent model is
empty. Models can then be built by using [model
modification](#model-modification) methods to add columns and rows.

### Model solution

#### Solving a model

HiGHS has a single method for solving the incumbent model. By default,
HiGHS will use what it considers to be the best solver and run-time
option settings for the incumbent model. To force HiGHS to use a
particular solver and/or run-time options, [`HighsOptions`](#options)
values can be modified. After this is run, the solution data that is
available from HiGHS comes from four sources: the model status, a set
of scalar information values, the solution of the model and the status
of each variable and constraint.

#### Model status

A single scalar `HighsModelStatus` value tells the user which of the following holds

* The solution is optimal
* The model is infeasible or unbounded
* A user-defined time or iteration limit has been reached
* A user-defined objective target or bound has been reached
* Some error has occurred

#### Scalar information

Scalar 'HighsInfo' values tell the user about

* The number of solver iterations performed
* Whether a primal or dual solution is available, and whether it is feasible
* Whether _basis_ information is valid
* The objective function value
* The number, maximum and sum of primal and dual infeasibilities
* The node count, dual bound and gap when solving MIP problems

#### Solution values

Arrays of primal and dual `HighsSolution` values are made available to
the user, if known. When solving MIP problem, no dual solution is
available, and if the model is found to be logically infeasible then
no solution values are available.

#### Variable and constraint status

Arrays of `HighsBasis` constraint status values are made available to
the user, if known. No constraint status values are available when
solving MIP problems, if the model is found to be logically
infeasible, or if the interior point solver is used without crossover.

### Model extraction

The incumbent model can be extracted from HiGHS.

### Model modification

Variables and constraints can be added to the incumbent model by
supplying the necessary costs, bounds and constraint matrix
entries. The Hessian matrix can be replaced. Costs, bounds and
constraint matrix entries can be changed. Variables and constraints
can be deleted from the incumbent model. When LP problems are
modified, the HiGHS simplex solvers will start from an advanced basis
whenever possible.

### Other features

### Options





