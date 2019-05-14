TreeMerge
=======
TreeMerge is a tool for scaling phylogengy estimation methods to large datasets. TreeMerge can be used in a divide-and-conquer framework as follows: 1) divide the species set into disjoint subsets, 2) construct trees on each subset, and 3) combine the subset trees using an associated distance matrix (on the full species set). TreeMerge has been sucessfully tested in the context of species tree estimation (Molloy and Warnow, 2019).


REQUIREMENTS
------------
+ Python 2.7 or later
+ [DendroPy](https://www.dendropy.org) 4.3.0
+ [PAUP*](http://phylosolutions.com/paup-test/) 4


TUTORIAL
--------
+ [Slides](http://erinkmolloy.web.illinois.edu/ekm-trees19.pdf)
+ [Hands-on Worksheet](https://github.com/ekmolloy/trees-in-the-desert-tutorial)


CITATION
--------
```
@article{MolloyWarnow2019,
    author={Erin Molloy and Tandy Warnow},
    title={{TreeMerge: A new method for improving the scalability of species tree estimation methods}},
    year={2019},
    journal={Bioinformatics},
    note={in press}
}
```

LICENSE
-------
Please see [LICENSE.txt](LICENSE.txt) for details.
