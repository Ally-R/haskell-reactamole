-- |
-- Module      :  Bio.Reactamole.ArrChoice
-- Copyright   :  (c) TBD
-- License     :  TBD
--
-- Maintainer  :  TBD
-- Stability   :  TBD
-- Portability :  TBD
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
  ( -- * Conditional SFs
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
entangle :: SF (Bool, (a, b)) (Either a b)
entangle = arrSp (\(PairS x (PairS y z)) -> CondS x y z)

--------------------------------------------------------------------------------

-- | Split the input between the two argument 'SF's, retagging and merging
-- their outputs.
--
-- Note that this is in general not a functor.
(+++) :: SF a b -> SF c d -> SF (Either a c) (Either b d)
f +++ g = arrSp (\(CondS x y z) -> PairS x (PairS y z))
  >>> (idSF *** (f *** g))
  >>> arrSp (\(PairS x (PairS y z)) -> CondS x y z)

-- | Fanin: split the input between the two argument arrows and merge their
-- outputs.
(|||) :: SF a c -> SF b c -> SF (Either a b) c
f ||| g = (f +++ g) >>> mergeSF

-- | Feed marked inputs through the argument arrow, passing the rest through
-- unchanged to the output.
left :: SF a b -> SF (Either a c) (Either b c)
left f = f +++ idSF

-- | A mirror image of left.
right :: SF a b -> SF (Either c a) (Either c b)
right f = idSF +++ f

--------------------------------------------------------------------------------

-- | Merge a conditional bool signal.
--
-- Implmentation uses the equivalence of the following two statements:
--   * @w = if x then y else z@
--   * @w = (x && y) || (!x && z)@
mergeBool :: SF (Either Bool Bool) Bool
mergeBool = arrSp (\(CondS x y z) ->
  PairS (PairS x y) (PairS x z))
  >>> (andSF *** (first notSF >>> andSF))
  >>> orSF

-- | Merge a conditional real signal.
--
-- If the conditional input signal is "@if x then y else z@", then the merged
-- real has solution \( w = x*y + (1-x)*z \) where the Bool \( x \) is
-- interpreted as a real-value that is @1@ when True and @0@ when False.
--
-- __NOTE__: Since Bools are only approximately 0- or 1-valued, the merged
--           real signal is also only approximately equal to @y@ and @z@.
mergeReal :: SF (Either Double Double) Double
mergeReal = arrSp (\(CondS (BoolS x x') y z) ->
  PairS (PairS y (RealS [Term 1 [x]]))
        (PairS z (RealS [Term 1 [x']])))
  >>> (multSF *** multSF)
  >>> addSF

-- | Merge a conditional pair signal.
mergePair :: SF (Either (a, b) (a, b)) (a, b)
mergePair = arrSp (\(CondS x (PairS ya yb) (PairS za zb)) ->
  PairS (CondS x ya za) (CondS x yb zb))
  >>> (mergeSF *** mergeSF)

-- | Merge a conditional tuple 3 signal.
mergeTup3 :: SF (Either (a, b, c) (a, b, c)) (a, b, c)
mergeTup3 = arrSp (\(CondS x (Tup3S y z w) (Tup3S y' z' w')) ->
  PairS (CondS x y y') (CondS x (PairS z w) (PairS z' w')))
  >>> (mergeSF *** mergePair)
  >>> arrSp (\(PairS x (PairS y z)) -> Tup3S x y z)

-- | Merge a conditional tuple 4 signal.
mergeTup4 :: SF (Either (a, b, c, d) (a, b, c, d)) (a, b, c, d)
mergeTup4 = arrSp (\(CondS x (Tup4S y z w u) (Tup4S y' z' w' u')) ->
  PairS (CondS x y y') (CondS x (Tup3S z w u) (Tup3S z' w' u')))
  >>> (mergeSF *** mergeTup3)
  >>> arrSp (\(PairS x (Tup3S y z w)) -> Tup4S x y z w)

-- | Merge a conditional tuple 5 signal.
mergeTup5 :: SF (Either (a, b, c, d, e) (a, b, c, d, e)) (a, b, c, d, e)
mergeTup5 = arrSp (\(CondS x (Tup5S y z w u v) (Tup5S y' z' w' u' v')) ->
  PairS (CondS x y y') (CondS x (Tup4S z w u v) (Tup4S z' w' u' v')))
  >>> (mergeSF *** mergeTup4)
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
mergeCond :: SF (Either (Either a b) (Either a b)) (Either a b)
mergeCond = arrSp (\(CondS x (CondS y y' y'') (CondS z z' z'')) ->
  PairS (CondS x y z) (PairS (CondS x y' z')
                             (CondS x y'' z'')))
  >>> (mergeBool *** (mergeSF *** mergeSF))
  >>> arrSp (\(PairS x (PairS y z)) -> CondS x y z)

-- | Merge a conditional signal.
mergeSF :: SF (Either c c) c
mergeSF = SF $ \s@(Sg (CondS x y z) sys ic) -> case (x,y,z) of
  (TrueS,  _, _) -> Sg y sys ic
  (FalseS, _, _) -> Sg z sys ic
  (_, NullS      , _) -> nullSg
  (_, TrueS      , _) -> runSF mergeBool s
  (_, FalseS     , _) -> runSF mergeBool s
  (_, BoolS _ _  , _) -> runSF mergeBool s
  (_, RealS _    , _) -> runSF mergeReal s
  (_, PairS _ _  , _) -> runSF mergePair s
  (_, Tup3S {}   , _) -> runSF mergeTup3 s
  (_, Tup4S {}   , _) -> runSF mergeTup4 s
  (_, Tup5S {}   , _) -> runSF mergeTup5 s
  (_, CondS {}   , _) -> runSF mergeCond s
