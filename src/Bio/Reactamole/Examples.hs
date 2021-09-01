-- |
-- Module      :  Bio.Reactamole.Examples
-- Copyright   :  (c) DigMP Research Group 2021
-- License     :  MIT
--
-- Examples of uses for Reactamole.
--
--  * Including all examples presented in:
--    /Reactamole: Functional Reactive Molecular Programming/, by Kinge,
--    Lathrop, Osera and Rogers in /DNA Computing and Molecular Programming/ 27,
--    September 2021.

module Bio.Reactamole.Examples
  ( -- * Example CRNs
    sinCRN
  , lowPass
  , bandPass
  , modulate
  , demodulate
  , srLatch
  , clock
  , unanimousCRN

    -- * Helpers
  , constMult
  , carrier
  , entangle
  , rectify
  , isPos
  , unanimous
  ) where

import Bio.Reactamole
import Bio.Reactamole.ArrChoice

-- | CRN for the sine function.
sinCRN :: CRN a Double
sinCRN = loop $ proj2 >>> negateCRN >>> integrate 1 >>> integrate 0 >>> dupCRN

-- | Multiply a 'Signal' by the given constant.
constMult :: Double -> CRN Double Double
constMult d = constRl d &&& idCRN >>> multCRN

-- lowPass1 :: CRN Double Double
-- lowPass1 = loop (second negateCRN >>> addCRN >>> intCRN >>> dupCRN)

-- | A lowpass filter.
lowPass :: Double -> Double -> CRN Double Double
lowPass a b = loop (constMult a *** constMult (-b) >>> addCRN >>> integrate 0 >>> dupCRN)


-- lowPass2 :: Double -> Double -> CRN Double Double
-- lowPass2 k w = loop (constMult k *** constMult (-w) >>> addCRN >>> intCRN >>> dupCRN)

-- bandPass1 :: CRN Double Double
-- bandPass1 = loop (second (negateCRN &&& intCRN >>> addCRN) >>> addCRN >>> intCRN >>> dupCRN)

-- bandPass2 :: Double -> Double -> Double -> CRN Double Double
-- bandPass2 k q w =
--   loop (second (negateCRN &&& intCRN >>> (constMult w *** constMult q >>> addCRN))
--         >>> first (constMult k) >>> addCRN >>> intCRN >>> dupCRN)

-- bandPass3 :: Double -> Double -> Double -> CRN Double Double
-- bandPass3 k q w = loop (first (constMult k)
--   >>> second (constMult (-w) &&& (intCRN >>> constMult q) >>> addCRN)
--   >>> addCRN >>> intCRN >>> dupCRN)

-- | A bandpass filter.
bandPass :: Double -> Double -> Double -> CRN Double Double
bandPass a b c = loop (first (constMult a)
  >>> second (constMult (-c) &&& (integrate 0 >>> constMult (-b)) >>> addCRN)
  >>> addCRN >>> integrate 0 >>> dupCRN)

-- | Generate a sine wave with the given frequency.
carrier :: Double -> CRN a Double
carrier w = loop (proj2 >>> constMult w >>> integrate 1
        >>> constMult (-w) >>> integrate 0
        >>> dupCRN)

-- | Produce a pair @m(t)@ that satisfies
-- \( \frac{dm}{dt} = u(t)\cdot\sin(ft) - m(t) \).
modulate :: Double -> CRN Double Double
modulate w = loop (first (idCRN &&& carrier w >>> multCRN)
  >>> second negateCRN >>> addCRN >>> integrate 0 >>> dupCRN)

-- | Demodulate a signal that has been modulated on a carrier signal at
-- frequency @w@.
demodulate :: Double -> Double -> CRN Double Double
demodulate w q = bandPass (w/q) (w/q) (w*w) >>> rectify >>> lowPass w w

-- expCRN :: CRN a Double
-- expCRN = loop (proj2 >>> intCRN >>> dupCRN)

-- tanh :: CRN a Double
-- tanh = loop (proj2 >>> square >>> negateCRN >>> plusOne >>> intCRN >>> dupCRN)
--   where square = dupCRN >>> multCRN
--         plusOne = constRl 1 &&& idCRN >>> addCRN

-- bandPass4 :: Double -> Double -> Double -> CRN Double Double
-- bandPass4 a b c =
--   loop
--     ( first (constMult a)
--         >>> second (constMult (- c) &&& (intCRN >>> constMult (- b)) >>> addCRN)
--         >>> addCRN
--         >>> intCRN
--         >>> dupCRN
--     )

rectify :: CRN Double Double
rectify = isPos &&& dupCRN >>> entangle >>> (idCRN ||| constRl 0)

-- | Determine if input double is positive.
isPos :: CRN Double Bool
isPos = posCRN

-- boolClock :: Double -> CRN () Bool
-- boolClock t = sin  >>> isPos

--------------------------------------------------------------------------------

-- | An SR latch CRN.
srLatch :: CRN (Bool, Bool) (Bool, Bool)
srLatch = loop $
  arrSp (\(PairS (PairS s' r') (PairS q q'))
         -> PairS (PairS s' q') (PairS r' q))
  >>> (nandCRN *** nandCRN)
  >>> dupCRN

-- | A clock CRN.
clock :: CRN a Bool
clock = sinCRN >>> isPos

-- | A Haskell function to determine if three inputs are the same.
--
-- For demonstration of lifting Boolean functions.
unanimous :: Bool -> Bool -> Bool -> Bool
unanimous x y z = if x then y && z else not (y || z)

-- | An CRN to determine if three inputs are the same.
--
-- Demonstrates the use of arr3Bl.
unanimousCRN :: CRN (Bool, Bool, Bool) Bool
unanimousCRN = arr3Bl unanimous
