---
title: Guide to HiGHS
permalink: /guide/
---

## Introduction

This guide gives abrief overview of what HiGHS can do when called from
an application. Detailed formal documentation setting out the data
structures and methods in the C++ interface, the corresponding methods
in the C interface, and how HiGHS can be controlled, is
[available](../docs/HighsDocumentation.pdf "HiGHS Documentation"). The
methods in the FORTRAN 90 and C# interfaces are functionally identical
to the methods in the C interface.

### Model defintion

Models are either passed to HiGHS directly, read from data files, or
built within HiGHS by adding columns and rows to an empty model. When
the model is passed or built, the matrix information is communicated
in [packed vector form](https://en.wikipedia.org/wiki/Sparse_matrix), either row-wise or column-wise. At any point in
time HiGHS has a single *incumbent model*.

##### Passing a model to HiGHS

When called from C++, a model is passed as an instance of the
`HighsModel` class (in general), or as an instance of the `HighsLp`
class if the model is a (mixed integer) linear programming problem. In
the C interface, the model is defined by passing the corresponding
scalar and array data as parameters.

##### Reading from a data file

HiGHS can read a model from files in the MPS or CPLEX LP formats. If
an MPS file contains a section corresponding to an extnesion not
supported by HiGHS, an error will be returned. For models in the CPLEX
LP format, HiGHS will read no further than the standard maximum line
length of 255 characters. HiGHS will return with an error if longer
lines are encountered.

##### Building a model

In a newly-created instance of HiGHS, the incumbent model is
empty. Models can then be built by using [model
modification](#model-modification) methods to add columns and rows.

### Model solution

##### Solving a model

HiGHS has a single method for solving the incumbent model. By default,
HiGHS will use what it considers to be the best solver and run-time
option settings for the incumbent model. To force HiGHS to use a
particular solver and/or run-time options, [`HighsOptions`
values](#options) may be modified. After this is run, the solution
data that is available from \HiGHS comes from four sources: the model
status, a set of scalar information values, the solution of the model
and the status of each variable and constraint.

##### Solving a model


### Model extraction

### Model modification

### Other features

### Options





