-- |
-- Module      :  Bio.Reactamole.ArrChoice
-- Copyright   :  (c) DigMP Research Group 2021
-- License     :  MIT
--
-- Arrow definitions for CRN conditionals and 'Bio.Reactamole.Core.Either'
-- values.
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

module Bio.Reactamole.ArrChoice
  ( -- * Conditional CRNs
    entangle
  , (+++)
  , (|||)
  , left
  , right
  ) where

import Bio.Reactamole.Bool
import Bio.Reactamole.Real
import Bio.Reactamole.Arr
import Bio.Reactamole.Core

--------------------------------------------------------------------------------

-- | The entangle signal function.
entangle :: CRN (Bool, (a, b)) (Either a b)
entangle = arrSp (\(PairS x (PairS y z)) -> CondS x y z)

--------------------------------------------------------------------------------

-- | Split the input between the two argument 'CRN's, retagging and merging
-- their outputs.
--
-- Note that this is in general not a functor.
(+++) :: CRN a b -> CRN c d -> CRN (Either a c) (Either b d)
f +++ g = arrSp (\(CondS x y z) -> PairS x (PairS y z))
  >>> (idCRN *** (f *** g))
  >>> arrSp (\(PairS x (PairS y z)) -> CondS x y z)

-- | Fanin: split the input between the two argument arrows and merge their
-- outputs.
(|||) :: CRN a c -> CRN b c -> CRN (Either a b) c
f ||| g = (f +++ g) >>> mergeCRN

-- | Feed marked inputs through the argument arrow, passing the rest through
-- unchanged to the output.
left :: CRN a b -> CRN (Either a c) (Either b c)
left f = f +++ idCRN

-- | A mirror image of left.
right :: CRN a b -> CRN (Either c a) (Either c b)
right f = idCRN +++ f

--------------------------------------------------------------------------------

-- | Merge a conditional bool signal.
--
-- Implmentation uses the equivalence of the following two statements:
--   * @w = if x then y else z@
--   * @w = (x && y) || (!x && z)@
mergeBool :: CRN (Either Bool Bool) Bool
mergeBool = arrSp (\(CondS x y z) ->
  PairS (PairS x y) (PairS x z))
  >>> (andCRN *** (first notCRN >>> andCRN))
  >>> orCRN

-- | Merge a conditional real signal.
--
-- If the conditional input signal is "@if x then y else z@", then the merged
-- real has solution \( w = x*y + (1-x)*z \) where the Bool \( x \) is
-- interpreted as a real-value that is @1@ when True and @0@ when False.
--
-- __NOTE__: Since Bools are only approximately 0- or 1-valued, the merged
--           real signal is also only approximately equal to @y@ and @z@.
mergeReal :: CRN (Either Double Double) Double
mergeReal = arrSp (\(CondS (BoolS x x') y z) ->
  PairS (PairS y (RealS [Term 1 [x]]))
        (PairS z (RealS [Term 1 [x']])))
  >>> (multCRN *** multCRN)
  >>> addCRN

-- | Merge a conditional pair signal.
mergePair :: CRN (Either (a, b) (a, b)) (a, b)
mergePair = arrSp (\(CondS x (PairS ya yb) (PairS za zb)) ->
  PairS (CondS x ya za) (CondS x yb zb))
  >>> (mergeCRN *** mergeCRN)

-- | Merge a conditional tuple 3 signal.
mergeTup3 :: CRN (Either (a, b, c) (a, b, c)) (a, b, c)
mergeTup3 = arrSp (\(CondS x (Tup3S y z w) (Tup3S y' z' w')) ->
  PairS (CondS x y y') (CondS x (PairS z w) (PairS z' w')))
  >>> (mergeCRN *** mergePair)
  >>> arrSp (\(PairS x (PairS y z)) -> Tup3S x y z)

-- | Merge a conditional tuple 4 signal.
mergeTup4 :: CRN (Either (a, b, c, d) (a, b, c, d)) (a, b, c, d)
mergeTup4 = arrSp (\(CondS x (Tup4S y z w u) (Tup4S y' z' w' u')) ->
  PairS (CondS x y y') (CondS x (Tup3S z w u) (Tup3S z' w' u')))
  >>> (mergeCRN *** mergeTup3)
  >>> arrSp (\(PairS x (Tup3S y z w)) -> Tup4S x y z w)

-- | Merge a conditional tuple 5 signal.
mergeTup5 :: CRN (Either (a, b, c, d, e) (a, b, c, d, e)) (a, b, c, d, e)
mergeTup5 = arrSp (\(CondS x (Tup5S y z w u v) (Tup5S y' z' w' u' v')) ->
  PairS (CondS x y y') (CondS x (Tup4S z w u v) (Tup4S z' w' u' v')))
  >>> (mergeCRN *** mergeTup4)
  >>> arrSp (\(PairS x (Tup4S y z w u)) -> Tup5S x y z w u)

-- | Merge a nested conditional signal.
--
-- Cosnider the following conditional:
--  * @if x then (if y then y' else y'')
--          else (if z then z' else z'')@
--
-- The above conditional is converted to: (@if w then w' else w''@) where
--  * @w@   is the merged Bool:      (@if x then y   else z  @)
--  * @w'@  is the merged Species a: (@if x then y'  else z' @)
--  * @w''@ is the merged Species b: (@if x then y'' else z''@)
mergeCond :: CRN (Either (Either a b) (Either a b)) (Either a b)
mergeCond = arrSp (\(CondS x (CondS y y' y'') (CondS z z' z'')) ->
  PairS (CondS x y z) (PairS (CondS x y' z')
                             (CondS x y'' z'')))
  >>> (mergeBool *** (mergeCRN *** mergeCRN))
  >>> arrSp (\(PairS x (PairS y z)) -> CondS x y z)

-- | Merge a conditional signal.
mergeCRN :: CRN (Either c c) c
mergeCRN = CRN $ \s@(Sg (CondS x y z) sys ic) -> case (x,y,z) of
  (TrueS,  _, _) -> Sg y sys ic
  (FalseS, _, _) -> Sg z sys ic
  (_, NullS      , _) -> nullSg
  (_, TrueS      , _) -> runCRN mergeBool s
  (_, FalseS     , _) -> runCRN mergeBool s
  (_, BoolS _ _  , _) -> runCRN mergeBool s
  (_, RealS _    , _) -> runCRN mergeReal s
  (_, PairS _ _  , _) -> runCRN mergePair s
  (_, Tup3S {}   , _) -> runCRN mergeTup3 s
  (_, Tup4S {}   , _) -> runCRN mergeTup4 s
  (_, Tup5S {}   , _) -> runCRN mergeTup5 s
  (_, CondS {}   , _) -> runCRN mergeCond s
