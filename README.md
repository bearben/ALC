ApproxLatCount (ALC) is a tool for lattices counting problem over linear constraints, i.e., in a convex polytope.
All experimental data including source codes of ALC and #SMT(LA)+ALC, and benchmarks, for AAAI2024 submission, can be found in "AAAI_Supplementary.7z".


ALC relies on Armadillo and GLPK. To compile ApproxLatCount, build on Ubuntu simply type:

```bash
sudo apt-get install g++
sudo apt-get install libglpk-dev
sudo apt-get install libboost-dev
sudo apt-get install libarmadillo-dev
make
```

And test ApproxLatCount:
> ./ApproxLatCount test/5_5_3.in <br />
> ./ApproxLatCount test/cube_10.in <br />
> ./ApproxLatCount test/simplex_8.in <br />
> ./ApproxLatCount test/timeo.in <br />

