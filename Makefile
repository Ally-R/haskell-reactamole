.PHONY : build
build :
	stack build

.PHONY : clean
clean :
	stack clean

.PHONY : repl
repl :
	stack repl

.PHONY : docs
docs :
	stack haddock --haddock-arguments "-o docs"
