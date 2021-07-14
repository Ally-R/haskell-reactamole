-- |
-- Module      :  Bio.Reactamole.Real
-- Copyright   :  (c) TBD
-- License     :  TBD
--
-- Maintainer  :  TBD
-- Stability   :  TBD
-- Portability :  TBD
--
-- Functions for manipulating Reals in Reactamole.

module Bio.Reactamole.Real
  ( -- * Signal Functions
    -- ** Basic SFs
    negateSF
  , constRl

    -- ** Comparators
  , posSF
  , ltSF
  , gtSF

    -- ** Real Function SFs
  , addSF
  , subSF
  , multSF
  , integrate
  ) where

import Bio.Reactamole.Arr
import Bio.Reactamole.Core
import Bio.Reactamole.Bool

--------------------------------------------------------------------------------

-- | Numerical negation.
negateSF :: SF Double Double
negateSF = arrSp (\(RealS eq) -> RealS (multEq eq [Term (-1) []]))

-- | Numerical addition for SFs.
addSF :: SF (Double, Double) Double
addSF = arrSp (\(PairS (RealS eq1) (RealS eq2)) -> RealS (addEq eq1 eq2))

-- | Numerical multiplication for SFs.
multSF :: SF (Double, Double) Double
multSF = arrSp (\(PairS (RealS eq1) (RealS eq2)) -> RealS (multEq eq1 eq2))

-- | Create an SF that ignores its input, producing the double d regardless.
constRl :: Double -> SF a Double
constRl d = SF $ const (Sg (RealS [Term d []]) [] [])

-- | Numerical subtraction for SFs.
subSF :: SF (Double, Double) Double
subSF = second negateSF >>> addSF

--------------------------------------------------------------------------------

-- | Numerical integration for SFs from zero to infinity.
integrate :: Double -> SF Double Double
integrate c = SF $ \(Sg (RealS eq) sys ic) ->
  let y  = length sys
      y' = y + 1
      posTerms = filter (\(Term c' _) -> c' > 0) eq
      negTerms = multEq [Term (-1) []] (filter (\(Term c' _) -> c' < 0) eq)
      eq' = [Term 1 [y], Term (-1) [y']]
      sys' = sys ++ [posTerms ++ [Term (-1) [y,y']], negTerms ++ [Term (-1) [y,y']]]
      y0 =  if c > 0 then c    else 0
      y'0 = if c < 0 then (-c) else 0
      ic' = ic ++ [y0, y'0]
  in Sg (RealS eq') sys' ic'

--------------------------------------------------------------------------------

-- | Check if a double is positive.
--
-- Gathers all the positive terms of the input and has them bias the Bool
-- towards True and gathers all the negative terms of the input and has them
-- bias the output towards False.
--
-- __WARNING__: Latency depends on how close @f(t)@ is to zero and therefore
--              may fail when @f(t)@ is close to zero.
posSF :: SF Double Bool
posSF = SF $ \(Sg (RealS eq) sys ic) ->
  let y  = length sys
      y' = y + 1
      posTerms = filter (\(Term c _) -> c > 0) eq
      negTerms = filter (\(Term c _) -> c < 0) eq
      yEq  = multEq [Term 100 [y']] posTerms ++ multEq [Term 100 [y]] negTerms
      y'Eq = multEq [Term (-1) []] yEq
  in Sg (BoolS y y') (sys ++ [yEq, y'Eq]) (ic ++ [0.5, 0.5])

-- | Check if @f(t) > g(t)@ where @f@ and @g@ are the two input signals.
--
-- __WARNING__: Latency depends on how close @f(t)@ and @g(t)@ are to one
--              another and therefore may fail when @f(t) - g(t)@ is small.
gtSF :: SF (Double, Double) Bool
gtSF = subSF >>> posSF

-- | Check if @f(t) < g(t)@ where @f@ and @g@ are the two input signals.
--
-- __WARNING__: Latency depends on how close @f(t)@ and @g(t)@ are to one
--              another and therefore may fail when @f(t) - g(t)@ is small.
ltSF :: SF (Double, Double) Bool
ltSF = gtSF >>> notSF
