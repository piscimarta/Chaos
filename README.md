## Building

Clone this repo with

```bash
git clone https://github.com/piscimarta/Chaos
```

or

```bash
git clone git@github.com:piscimarta/Chaos.git
```

### g++

You can compile and link with `g++`. For example, to link and compile the `main2.cpp` file:

```bash
g++ main2.cpp header.cpp system.cpp planet.cpp -I include -larmadillo -o main2
```
You might need to add the  `-std=gnu++11` if you are a Mac user.
You can then run the executable with

```bash
./main
```

