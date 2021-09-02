-- |
-- Module      :  Bio.Reactamole.Real
-- Copyright   :  (c) DigMP Research Group 2021
-- License     :  MIT
--
-- Functions for manipulating Reals in Reactamole.

module Bio.Reactamole.Real
  ( -- * Signal Functions
    -- ** Basic CRNs
    negateCRN
  , constRl

    -- ** Comparators
  , posCRN
  , ltCRN
  , gtCRN

    -- ** Real Function CRNs
  , addCRN
  , subCRN
  , multCRN
  , integrate
  ) where

import Bio.Reactamole.Arr
import Bio.Reactamole.Core
import Bio.Reactamole.Bool

--------------------------------------------------------------------------------

-- | Numerical negation.
negateCRN :: CRN Double Double
negateCRN = arrSp (\(RealS eq) -> RealS (multEq eq [Term (-1) []]))

-- | Numerical addition for CRNs.
addCRN :: CRN (Double, Double) Double
addCRN = arrSp (\(PairS (RealS eq1) (RealS eq2)) -> RealS (addEq eq1 eq2))

-- | Numerical multiplication for CRNs.
multCRN :: CRN (Double, Double) Double
multCRN = arrSp (\(PairS (RealS eq1) (RealS eq2)) -> RealS (multEq eq1 eq2))

-- | Create an CRN that ignores its input, producing the double d regardless.
constRl :: Double -> CRN a Double
constRl d = CRN $ const (Sg (RealS [Term d []]) [] [])

-- | Numerical subtraction for CRNs.
subCRN :: CRN (Double, Double) Double
subCRN = second negateCRN >>> addCRN

--------------------------------------------------------------------------------

-- | Numerical integration for CRNs from zero to infinity.
integrate :: Double -> CRN Double Double
integrate c = CRN $ \(Sg (RealS eq) sys ic) ->
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
posCRN :: CRN Double Bool
posCRN = CRN $ \(Sg (RealS eq) sys ic) ->
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
gtCRN :: CRN (Double, Double) Bool
gtCRN = subCRN >>> posCRN

-- | Check if @f(t) < g(t)@ where @f@ and @g@ are the two input signals.
--
-- __WARNING__: Latency depends on how close @f(t)@ and @g(t)@ are to one
--              another and therefore may fail when @f(t) - g(t)@ is small.
ltCRN :: CRN (Double, Double) Bool
ltCRN = gtCRN >>> notCRN
