## Most Parsimonious Likelihood (has lumps)

You need to have `sage` installed for this, which you can get [here](http://www.sagemath.org/download-linux.html), as well as using
```
sudo apt-get install sagemath
```
Then you run it by running `sage` in interactive mode and calling `load("mpl.sage")`, which will print
```
[(v_1, 0), (v_2, 0), (v_3, 0), (v_4, 0)] [(v_u, 0), (v_v, 0)]
[(v_1, 1), (v_2, 0), (v_3, 0), (v_4, 0)] [(v_u, 0), (v_v, 0)]
[(v_1, 0), (v_2, 1), (v_3, 0), (v_4, 0)] [(v_u, 0), (v_v, 0)]
[(v_1, 0), (v_2, 0), (v_3, 1), (v_4, 0)] [(v_u, 0), (v_v, 0)]
[(v_1, 0), (v_2, 0), (v_3, 0), (v_4, 1)] [(v_u, 0), (v_v, 0)]
[(v_1, 1), (v_2, 1), (v_3, 0), (v_4, 0)] [(v_u, 1), (v_v, 0)]
The usual example we use is in the variable data = [3, 1, 1, 1, 1, 2]
Here's some random data: random_data = [3, 4, 3, 6, 7, 4]

Usage:
    summarize_lumps(data)
    int_max, border_max = sort_functions(data); find_lumps(int_max, border_max)
    
The sort_functions procedure takes a long time, so you can get intermediate values.

    print_lumps(find_lumps(int_max, border_max))  # prints LaTeX tables. 

You can make your own data, it just needs to be a list of length 6.
```
You can also switch to a different tree by typing 
```
g = tree5
```
or 
```
g = tree3
```
but you need to make your own if you want some other trees. Printing out ```g``` will provide helpful information for this.
