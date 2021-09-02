# Reactamole: A Functional Reactive Molecular Programming DSL

![CI Workflow Badge](https://github.com/digMP/haskell-reactamole/actions/workflows/main.yml/badge.svg)

Reactamole is a domain-specific language for molecular programming that
utilizes reactive functional programming principles. Reactamole observes a
direct correspondence between signal functions in a reactive functional program
(FRP) and chemical reaction networks (CRNs) in a molecular program. This
correspondence allows us to directly translate the core combinators of
(arrowized) FRP to CRNs. Because Reactamole is embedded in the Haskell
programming language, takes advantage of Haskell's rich type system to
ensure the well-formedness of the resulting chemical reaction networks.

## Installation

Reactamole is implemented as a [Haskell Stack](https://haskellstack.org)
project and as such requires a Haskell toolchain (namely `stack`) to build and
run.

1.  Download `ghcup` (https://www.haskell.org/ghcup) to install Haskell and its
    toolchain.
2.  `stack build` to build the project.
3.  `stack repl` to run Reactamole in GHCi.
4.  `stack haddock --haddock-arguments "-o docs"` to build API documentation
    (deposited in `/docs`).

The included `Makefile` also includes shortcuts for these commands if you are
unfamiliar with working within the Haskell ecosystem.

## Running Examples

The `Bio.Reactamole.Examples` module contains a number of example molecular
programs for you to explore. You can use the functions from
`Bio.Reactamole.Export` to export these Reactamole programs as ODEs or
collections of reactions. Here is an example of their usage:

~~~console
$> stack repl
...
GHCi, version 8.8.4: https://www.haskell.org/ghc/  :? for help
...
Ok, 10 modules loaded.
Loaded GHCi configuration from /private/var/folders/rj/4hpzks9x6m3bvbdbvqbs3ndm0000gn/T/haskell-stack-ghci/38c8115e/ghci-script
λ> import Bio.Reactamole
λ> import Bio.Reactamole.Examples
λ> :t srLatch
srLatch :: CRN (Bool, Bool) (Bool, Bool)
λ> toIVP srLatch
INPUT:
  Tuple(Bool(x0, x1), Bool(x2, x3))

OUTPUT:
  Tuple(Bool(x4, x5), Bool(x6, x7))

EQUATIONS:
  dx4/dt = [-1.0e-2*x0*x4*x6,+1.0e-2*x1*x5,+3.0e-2*x4^2*x5,-3.0e-2*x4*x5^2,+1.0e-2*x5*x7]
  dx5/dt = [+1.0e-2*x0*x4*x6,-1.0e-2*x1*x5,-3.0e-2*x4^2*x5,+3.0e-2*x4*x5^2,-1.0e-2*x5*x7]
  dx6/dt = [-1.0e-2*x2*x4*x6,+1.0e-2*x3*x7,+1.0e-2*x5*x7,+3.0e-2*x6^2*x7,-3.0e-2*x6*x7^2]
  dx7/dt = [+1.0e-2*x2*x4*x6,-1.0e-2*x3*x7,-1.0e-2*x5*x7,-3.0e-2*x6^2*x7,+3.0e-2*x6*x7^2]

INITIAL CONDITIONS:
  x4(0) = 0.0
  x5(0) = 1.0
  x6(0) = 0.0
  x7(0) = 1.0

λ> toRxns srLatch
INPUT:
  Tuple(Bool(x0, x1), Bool(x2, x3))

OUTPUT:
  Tuple(Bool(x4, x5), Bool(x6, x7))

REACTIONS:
  x0 + x4 + x6 --{1.0e-2}-> x0 + x6 + x5
  x1 + x5 --{1.0e-2}-> x1 + x4
  x2 + x4 + x6 --{1.0e-2}-> x2 + x4 + x7
  x3 + x7 --{1.0e-2}-> x3 + x6
  x4 + x4 + x5 --{3.0e-2}-> x4 + x4 + x4
  x4 + x5 + x5 --{3.0e-2}-> x5 + x5 + x5
  x5 + x7 --{1.0e-2}-> x4 + x6
  x6 + x6 + x7 --{3.0e-2}-> x6 + x6 + x6
  x6 + x7 + x7 --{3.0e-2}-> x7 + x7 + x7

INITIAL CONDITIONS:
  x4(0) = 0.0
  x5(0) = 1.0
  x6(0) = 0.0
  x7(0) = 1.0
~~~

## Publication

Titus H. Klinge, James I. Lathrop, Peter-Michael Osera, and Allison Rogers.
Reactamole: Functional Reactive Molecular Programming. 27th International
Conference on DNA Computing and Molecular Programming (DNA '27). September,
2020, Oxford, UK. DOI:10.4230/LIPIcs.CVIT.2016.23.
