### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ aff9313b-2fd8-46ca-b399-47fa21b324c5
begin
    import Pkg
    Pkg.activate(Base.current_project())
    Pkg.instantiate()
end

# ╔═╡ bcca4c48-89be-11ec-3c89-3189fd524f08
md"""
# Template-matching in Julia using Pluto.jl

## Introduction

Julia is a general programming language suited for scientific computing applications. Moreover, Julia is a dynamic, high-level language: executing a program written in Julia does not require the user to compile it beforehand, as for C or Fortran, and its syntax is similar to MATLAB or Python. Nonetheless, its performance is close to ones of static languages.

Pluto is a Julia library for writing notebooks similar (but better in many respect) to Jupyter notebooks. In this notebook I will work out an example and very briefly introduce Julia and Pluto. The goal is to template match some signals from acoustic measurements of hydraulic fracture growth performed under triaxial confinement in laboratory (more informations [here](https://www.researchgate.net/profile/Seyyedmaalek-Momeni/publication/352780078_Combining_active_and_passive_acoustic_methods_to_image_hydraulic_fracture_growth_in_laboratory_experiments/links/60d87837a6fdccb745ea29dd/Combining-active-and-passive-acoustic-methods-to-image-hydraulic-fracture-growth-in-laboratory-experiments.pdf)). 
"""

# ╔═╡ 11ba1299-05e6-46c1-8716-ac06d5fc46a5
md"""
## Julia and Pluto.jl

### Pluto.jl

A Pluto notebook is a sequence of cells. Each cell contain some Julia code, which can be executed. In Julia, every statement is an expression and hence evaluates to a value. Pluto converts that value to an HTML element and shows it. What you are reading is the output of a cell containing just a string, which returns the string itself [^1]. Pluto is able to parse the content of cells and detect dependencies between one cell and another. When the user evaluates a cell, for example changing the value of a variable, Pluto will re-evaluate all the cells that depend on that one and in the right order.


Notice that Pluto requires one statement per cell, hence if you want to put more than one statement in a single cell you need to wrap them between the keywords `begin` and `end`, which are used to chain statements. The resulting block will evaluate to the value of the last statement.

For example

```julia
begin
	x = 1
	y = x + 1
	"hello world"
	2 * y
end
```

evaluates to 4.

[^1]: Actually, the value that is rendered as the piece of text you are reading is not a string, but an object whose type is defined in the Markdown standard library. You can create these object prefixing the string literals with `md`. In this way, the content will be formatted using Markdown. Otherwise you can use plain strings.
"""

# ╔═╡ 6dca0719-eeb3-42d5-8585-f89a77e579b3
md"""
### Defining Functions in Julia

Julia does not have methods attached to objects and accessible using the dot notation like in Python or C++. The main way of organizing code in Julia is by writing functions, and eventually grouping them into modules.

Function are defined using the `function` keyword. For example, the following will multiply the first two arguments and add the product to the third.

```julia
function f(x, y, z)
	x * y + z
end
```

Th `return` keyword is implied in the last statement of a function definition, but can be supplyied if one wants to.

Notice that the first time a function is called, Julia compiles on the fly that function. This is called just in time compilation or JIT. Hence the execution time of a function call change drastically between the first and subsequent calls. This is something to keep in mind when writing Julia code. There are some subtleties related to JIT compilation and the types of the arguments of a function, but we will talk about that in the next section.
"""

# ╔═╡ 364a442a-17ef-455d-9608-aecd3126db59
md"""
### Arrays and other types in Julia

In theory, one might write Julia code without thinking or caring about types, like in Python. However, it is practically impossibile to write and debug efficient code without some knowledge of the topic. One of the most important things to know is that the compiler can produce optimized machine code only for type-stable functions, where a function is called type-stabile if the type of its output is determined only from the types of its arguments (and not from their values). 

The most used type in scientific computing is probably the array, i.e. a data type that represent indexed objects of the same type and contiguous in memory. Julia has a rich type system that includes arrays. In particular, arrays are an example of parametric types in Julia, i.e. types that depends on other types. The arguments of a dependent type are surrounded by braces, so for example `Array{Float64, 1}` is the type of one-dimensional arrays of 64-bit floating points numbers. To produce an exemplar of this type, you have to call its constructor. The constructor is a function with the same name of the type. In the case of one-dimensional array (or vector) constructors, they take two arguments, the value with which filling the array and the number of elements, so for example `Array{Float64, 1}(1.0, 10)` will be a ten element vector filled with the number 1.0. Array indexing in Julia is similar to Python with one important exceptions, in Julia one starts to count from 1 instead of 0, as you would do in Fortran or Matlab. Also, the last element is indicated with the keyword `end` and you can apply a function that acts on scalars element-wise on arrays posfixing a dot to its name in the function call, as in the function definition below (really, the dot notation is used for broadcasting, a powerful concept that should be familiar to Numpy user).
"""

# ╔═╡ Cell order:
# ╠═aff9313b-2fd8-46ca-b399-47fa21b324c5
# ╟─bcca4c48-89be-11ec-3c89-3189fd524f08
# ╟─11ba1299-05e6-46c1-8716-ac06d5fc46a5
# ╟─6dca0719-eeb3-42d5-8585-f89a77e579b3
# ╟─364a442a-17ef-455d-9608-aecd3126db59
