-- |
-- Module      :  Bio.Reactamole.Arr
-- Copyright   :  (c) TBD
-- License     :  TBD
--
-- Maintainer  :  TBD
-- Stability   :  TBD
-- Portability :  TBD
--
-- Arrow definitions for CRNs in Reactamole.
--
-- Sources:
--
--  * /Generalising Monads to Arrows/, by John Hughes,
--    /Science of Computer Programming/ 37, pp67-111, May 2000.
--
--  * /A New Notation for Arrows/, by Ross Paterson, in /ICFP 2001/,
--    Firenze, Italy, pp229-240.
--
--  * [Control.Arrow package.](https://hackage.haskell.org/package/base-4.15.0.0/docs/Control-Arrow.html)

module Bio.Reactamole.Arr
  ( -- * Arrows for CRNs
    (<<<)
  , (>>>)
  , (***)
  , (&&&)
  , first
  , second
  , loop
  --, loop'

    -- * Species Instantiation
  , instSF

  --, derivativeOfTerm
  --, derivativeOfEq
  --, dropSg
  ) where

import Bio.Reactamole.Core

--------------------------------------------------------------------------------

-- | Right-to-left composition.
(<<<) :: SF b c -> SF a b -> SF a c
(SF f) <<< (SF g) = SF (f . g)

-- | Left-to-right composition.
(>>>) :: SF a b -> SF b c -> SF a c
f >>> g = g <<< f

-- | Split the input between the two argument arrows and combine their output.
--
-- Note that this is in general not a functor.
(***) :: SF a b -> SF a' b' -> SF (a, a') (b, b')
(SF f) *** (SF g) = SF $ \s@(Sg _ sys ic) ->
  let n = length sys
      (Sg spf sysf icf) = (dropSg n . f . runSF fstSF) s
      k = length sysf
      (Sg spg sysg icg) = (shiftSg n k . dropSg n . g . runSF sndSF) s
    in Sg (PairS spf spg) (sys ++ sysf ++ sysg) (ic ++ icf ++ icg)

-- | Fanout: send the input to both argument arrows and combine their output.
(&&&) :: SF a b -> SF a c -> SF a (b, c)
f &&& g = dupSF >>> (f *** g)

-- | Send the first component of the input through the argument arrow, and copy
-- the rest unchanged to the output.
first :: SF a b -> SF (a, c) (b, c)
first f = f *** idSF

-- | A mirror image of 'first'.
second :: SF a b -> SF (c, a) (c, b)
second f = idSF *** f

--------------------------------------------------------------------------------

-- | Loop the second parameter, c, back into itself.
--
-- The c in the resulting SF is unused and unmodified by the new SF
-- (essentially the identity).
loop' :: SF (a,c) (b,c) -> SF (a,c) (b,c)
loop' f = SF $ \s@(Sg (PairS _ y) _ _) ->
  let s'@(Sg (PairS _ y') _ _) = runSF (second instSF) s
      s''@(Sg (PairS _ y'') _ _) =  runSF (f >>> second instSF) s'
      Sg (PairS x _) sys ic = reWireVars s'' $ zip (varsSp y'') (varsSp y')
  in Sg (PairS x y) sys ic

-- | Feed the second component back in on itself.
loop :: HasDefault c => SF (a, c) (b, c) -> SF a b
loop f = SF $ \(Sg x sys ic) ->
  let SF g = loop' f
      Sg (PairS y _) sys' ic' = g (Sg (PairS x getDefault) sys ic)
  in Sg y sys' ic'

--------------------------------------------------------------------------------

-- | Instantiate all /uninstantiated/ 'Species'.
instSF :: SF a a
instSF = SF $ \s@(Sg sp _ _) -> case sp of
  NullS     -> s
  TrueS     -> runSF instBl s
  FalseS    -> runSF instBl s
  BoolS _ _ -> runSF instBl s
  RealS _   -> runSF instRl s
  PairS _ _ -> runSF (instSF *** instSF) s
  Tup3S {}  -> runSF instTup3 s
  Tup4S {}  -> runSF instTup4 s
  Tup5S {}  -> runSF instTup5 s
  CondS {}  -> runSF instCond s

-- | Instantiate 'BoolS' species.
instBl :: SF Bool Bool
instBl = SF $ \(Sg sp sys ic) ->
  let n = length sys
      constEq = [Term 0.0 []]
  in case sp of
    TrueS     -> Sg (BoolS n (n+1)) (sys ++ [constEq, constEq]) (ic ++ [1,0])
    FalseS    -> Sg (BoolS n (n+1)) (sys ++ [constEq, constEq]) (ic ++ [0,1])
    BoolS _ _ -> Sg sp sys ic

-- | Return derivative of the 'Term'.
derivativeOfTerm :: [Equation] -> Term -> Equation
derivativeOfTerm _ (Term _ []) = []
derivativeOfTerm sys (Term c (x:xs)) =
  let x'  = sys !! x
      xs' = derivativeOfTerm sys (Term 1 xs)
  in addEq (multEq x' [Term c xs]) (multEq [Term c [x]] xs')

-- | Return derivative of the 'Equation'.
derivativeOfEq :: [Equation] -> Equation -> Equation
derivativeOfEq sys = concatMap (derivativeOfTerm sys)

-- | Instantiate 'RealS' species.
instRl :: SF Double Double
instRl = SF $ \(Sg (RealS eq) sys ic) ->
  case eq of
    [Term 1 [_], Term (-1) [_]] -> Sg (RealS eq) sys ic
    [Term (-1) [_], Term 1 [_]] -> Sg (RealS (reverse eq)) sys ic
    _ -> let x     = length sys
             x'    = x + 1
             eq'   = [Term 1 [x], Term (-1) [x']]
             eqdt  = derivativeOfEq sys eq
             dxdt  = filter (\(Term c _) -> c > 0) eqdt
             dx'dt = filter (\(Term c _) -> c < 0) eqdt
             sys'  = if null dxdt && null dx'dt
                      then sys ++ [dxdt, dx'dt]
                      else sys ++ [dxdt  ++ [Term (-1) [x, x']],
                                   dx'dt ++ [Term (-1) [x, x']]]
             res   = evalEq eq ic
             x0    = if res > 0 then res  else 0
             x'0   = if res < 0 then -res else 0
             ic'   = ic ++ [x0, x'0]
          in Sg (RealS eq') sys' ic'

-- | Instantiate 'Tup3S' species.
instTup3 :: SF (a,b,c) (a,b,c)
instTup3 = arrSp (\(Tup3S x y z) -> PairS x (PairS y z))
             >>> (instSF *** instSF)
             >>> arrSp (\(PairS x (PairS y z)) -> Tup3S x y z)

-- | Instantiate 'Tup4S' species.
instTup4 :: SF (a,b,c,d) (a,b,c,d)
instTup4 = arrSp (\(Tup4S x y z w) -> PairS x (Tup3S y z w))
             >>> (instSF *** instSF)
             >>> arrSp (\(PairS x (Tup3S y z w)) -> Tup4S x y z w)

-- | Instantiate 'Tup5S' species.
instTup5 :: SF (a,b,c,d,e) (a,b,c,d,e)
instTup5 = arrSp (\(Tup5S x y z w u) -> PairS x (Tup4S y z w u))
             >>> (instSF *** instSF)
             >>> arrSp (\(PairS x (Tup4S y z w u)) -> Tup5S x y z w u)

-- | Instantiate 'CondS' Species.
instCond :: SF (Either a b) (Either a b)
instCond = arrSp (\(CondS x y z) -> PairS x (PairS y z))
             >>> (instSF *** instSF)
             >>> arrSp (\(PairS x (PairS y z)) -> CondS x y z)

--------------------------------------------------------------------------------

-- | Drop the first n variables of a signal.
--
-- Remove the first n _equations_ and the first n _initial conditions_ from the
-- initial value problem
--
-- __WARNING__: This does not re-index the variables! Use with caution!
--              Consider using with 'shiftSg' to re-index them.
dropSg :: Int -> Signal a -> Signal a
dropSg n (Sg sp sys ic) = Sg sp (drop n sys) (drop n ic)
