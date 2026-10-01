## Important note 

**This is a new implementation of ffiasm. If you need access to the original (and now obsoleted) implemenation at  https://github.com/iden3/ffiasm-old .**

# ffiasm

This package is a script that generates a Finite field Library in Intel64 / ARM64 Assembly

## Usage

install g++ nasm ang gmp library if you don't have it.

```
npm install -g ffiasm
mkdir myProject
cd myProject
buildzqfield -q 21888242871839275222246405745257275088548364400416034343698204186575808495617 -n Fr
```

You now will have the files fr.asm, fr.hpp, fr.cpp, fr_element.hpp, fr_generic.cpp, fr_raw_generic.cpp and fr_raw_arm64.s

```
ls
```

If you are in an Intel64 machine, assemble fr.asm. In linux:

```
nasm -felf64 fr.asm
```

In a mac:

```
nasm -fmacho64 --prefix _ fr.asm
```

Create a file named main.cpp to use the library

```C
#include <stdio.h>
#include <stdlib.h>
#include "fr.hpp"

int main() {
    Fr_init();

    FrElement a;
    a.type = Fr_SHORT;
    a.shortVal = 2;

    FrElement b;
    b.type = Fr_SHORT;
    b.shortVal = 6;

    FrElement c;

    Fr_mul(&c, &a, &b);

    char *c1 = Fr_element2str(&c);
    printf("Result: %s\n", c1);
    free(c1);
}
```

Compile it. In Intel64:

```
g++ -DUSE_ASM -DARCH_X86_64 main.cpp fr.o fr.cpp -o example -lgmp
```

In ARM64:

```
g++ -DUSE_ASM -DARCH_ARM64 main.cpp fr.cpp fr_generic.cpp fr_raw_generic.cpp fr_raw_arm64.s -o example -lgmp
```

Or, in any architecture, without assembly:

```
g++ main.cpp fr.cpp fr_generic.cpp fr_raw_generic.cpp -o example -lgmp
```

Run it
```
./example
```

# Benchmark

```
npm run benchmark
```

## License

ffiasm is part of the iden3 project copyright 2020 0KIMS association and published with LGPL-3 license. Please check the COPYING file for more details.

