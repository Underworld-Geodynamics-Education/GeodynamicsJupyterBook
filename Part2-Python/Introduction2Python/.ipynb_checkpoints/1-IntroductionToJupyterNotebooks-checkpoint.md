---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: '0.9'
    jupytext_version: 1.5.2
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

+++ {"deletable": true, "editable": true}

# Jupyter Notebooks 

We are going to be using `jupyter` notebooks in our examples a version of python known as `ipython` in the notebooks. 

The notebook is a __live__ document that can contain a number of different types of content, 
including fragments of programs that can actually be run within the text. It requires some form
of computational substrate to run (otherwise you will just see the preview). 

You can use the notebooks as a scratchpad to process data, for example, 
and the workflow that you develop is recorded and can be revisited. 
It's like a lab notebook but for data manipulation and programming. 
You can export the document you produce as a pdf to make a lab report (etc).

+++

## Notebook cells

A notebook is made up of different __cells__ which contain different types of content and which can be _executed_ by selecting the cell and _running_ it (see toolbar).
The results of running a cell depend on the content of the cell. It may produce an output cell beneath it with results, or it may reformat itself.

Try double-clicking on this cell and you will see the raw programming language underneath that produces the formatted text. 
If you run this cell it will display the text again.
The cells can have different types of content and will behave accordingly when run. 

````{note}

### Formatted text

This uses an almost-already-formatted version of text known as `markdown` to make content that can be formatted to look similar to a word document, but which also tends to highlight the intended formatting even if not processed. It is quite a useful form to learn for taking notes.

If you look at the raw content you will see how to make

   - bullet points
     - nested bullet points as well as 
     - text in __bold__ 
     - text in _italics_ (or emphasised) 
   
``` sh
echo "Code that can be highlighted for different languages"
ls -l
```

```python
print ("including python, unix shell scripts, fortran")

for i in range(0,100):
    print i
```

````


### Formatted text

This uses an almost-already-formatted version of text known as `markdown` to make content that can be formatted to look similar to a word document, but which also tends to highlight the intended formatting even if not processed. It is quite a useful form to learn for taking notes.

If you look at the raw content you will see how to make

   - bullet points
     - nested bullet points as well as 
     - text in __bold__ 
     - text in _italics_ (or emphasised) 
   
``` sh
echo "Code that can be highlighted for different languages"
ls -l
``````

```python
print ("including python, unix shell scripts, fortran")

for i in range(0,100):
    print i
```

Mathematical symbols and formulae that can appear in the text like this _"the circumference of a circle is $2\pi r$"_ or as equations like this:

$$
\begin{equation}
    A = \pi r^2
\end{equation}
$$

BUT you will need to be able to parse the $\LaTeX$ language to write mathematics. 

You can add links like this: [also see the markdown cheatsheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet)

And images like this: 

Running cells containing formatted code will replace the cell with the formatted version. This allows you to write a well-formatted page with interleaved code which also can be executed.


### Python code

Cells which are defined to be code cells contain python statements which will be executed when the cell is run. If multiple cells have code, then they will be run in the order you choose (or if you "run all" or "run all above" they will be run in the order they are listed).

```{code-cell} ipython3
---
jupyter:
  outputs_hidden: false
solution: hidden
---
## Example of a simple python code cell

print( "Hello little world")
a = 1

## The last statement in a cell prints its value 
a

## (this is sometimes a little confusing - add a pass statement to get rid of this !)

#pass
```

+++ {"deletable": true, "editable": true}

When you run a code cell, it is just the same as typing all the code into the interpreter. If you run a cell twice it is the same as if you typed everything twice ... the notebook remembers ! You can run some of a notebook or all of it, and you can run cells out of order. This is great for experimenting and getting things right, but be careful, this can break things easily.

```{code-cell} ipython3
---
jupyter:
  outputs_hidden: false
---
print ("Run number {}".format(a))
        
a += 1
```

+++ {"deletable": true, "editable": true}

### Hidden cells

In this notebook, I have added some extensions that mean some cells are hidden (I use this to hide some of the details when we learn something complicated). The hidden cells are small, blank lines that you can view by selecting them and hitting the carat button in the toolbar. 

There is one hidden below ... see if you can find it !

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
---
## The simplest possible python program

print ("You can run but you can't hide")
```

```{code-cell} ipython3
---
deletable: true
editable: true
hide_input: true
jupyter:
  outputs_hidden: false
solution: hidden
---
## This is a hidden cell !

## You can't usually see it but it still runs if you execute the notebook

print ("Yes you can !")
```

+++ {"deletable": true, "editable": true, "solution": "hidden", "solution_first": true}

### Exercise cells

Some cells are intended to be used for class __exercises__ and have a little plus sign which indicates you can expand the cell to see a worked example or answer to a question. These cells are not executed unless they are expanded and run individually. Try this one both ways.

```{code-cell} ipython3
---
deletable: true
editable: true
jupyter:
  outputs_hidden: false
solution: hidden
---
print ("You can hide and you can't run !")
```

+++ {"solution": "hidden", "solution_first": true}

## A more interesting exercise

Why does this not work ?

```python
print ("Run number {}".format(b) )
b += 1
```

and, more importantly, how would you fix this ?

```{code-cell} ipython3
---
jupyter:
  outputs_hidden: false
solution: hidden
---
# Try it !!
```

```{code-cell} ipython3
:solution: hidden

## Because b hasn't been defined. 

try:
    print ("Run number {}".format(c))
except: 
    print ("Run number 1")
    c = 1
        
c += 1
```