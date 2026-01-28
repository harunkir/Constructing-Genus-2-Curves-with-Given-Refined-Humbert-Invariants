# Constructing Genus 2 Curves with Given Refined Humbert Invariants

This repository contains the SageMath implementation of the algorithms described in the paper **"Constructing Genus 2 Curves with Given Refined Humbert Invariants"** by **Harun KIR**.

The code implements **Algorithm 23** from the paper, which constructs a principally polarized abelian surface $(A, \Theta)$ associated with a given geometric ternary quadratic form (in the sense of the paper).

## Requirements

* **SageMath**: The code was developed and tested using **SageMath 10.6**.

## Files

* **`(Simplified) Main algorithm`**: The main algorithm file. This script is self-contained and requires no external dependencies other than SageMath. It implements the code for both primitive and imprimitive geometric forms, and includes the verification examples given in the last section of the paper.

* **`Geometric form`**: This is an extra part of our codes. This code checks whether a given integral ternary quadratic form is a geometric form or not (in the sense of the paper). This is an extra because the termination of this code is not discussed in the paper. We just added this code here because one can still use this code if the code stops for the given (heuristic) bound, then it returns whether the given form is geometric or not.

## Usage

You can run the code by loading it into a SageMath session or executing it directly.

### Input Format

The algorithms utilize the standard SageMath convention for integral ternary quadratic forms. A list or tuple of 6 (integer) entries 

`[a, b, c, r, s, t]` corresponds to the quadratic form:

$$
q(x,y,z) = ax^2 + by^2 + cz^2 + ryz + sxz + txy
$$

The input must be a **geometric form** (a positive definite integral ternary quadratic form with specific properties defined in the paper).

### Running the Examples

The code includes the specific test cases discussed in the final section of the paper. These are:

* `[4, 12, 28, 0, 4, 4]` (Imprimitive)
* `[4, 4, 5, 4, 4, 4]` (Primitive)
* `[4, 4, 9, 4, 4, 4]` (Primitive)
* `[4, 4, 5, 0, 0, -4]` (Primitive)
* `[4, 4, 9, 0, 0, -4]` (Primitive)

Running the main algorithm on these inputs produces output similar to:

Given Ternary Form: [4, 12, 28, 0, 4, 4] [Imprimitive] D=-296, \tilde{q}=[11, 10, 9], s=(2, 446, 9) and verified True

Given Ternary Form: [4, 4, 5, 4, 4, 4] [Primitive] D=-11, \tilde{q}=[3, 1, 1], s=(2, 122, 9) and verified True.


The output provides the necessary construction parameters $(A, \Theta)$ as described in **Algorithm 23**. The `verified True` confirms that the constructed object is correct by using Proposition 6 of the  paper.

### Extensive Testing

In addition to the examples above, the file includes a list of 199 ternary quadratic forms from the author's [PhD thesis](https://queensu.scholaris.ca/server/api/core/bitstreams/b08499c4-5a74-4bb4-ae2b-c2f9d12b4d3a/content) (Chapter 8). The main algorithm successfully constructs and verifies $(A, \Theta)$ for all these cases. We also utilized the geometric form checker to verify that all these forms are indeed geometric.


### Licence

This project is licensed under the MIT License - see the LICENSE file for details.


