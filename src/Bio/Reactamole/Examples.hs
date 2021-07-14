-- |
-- Module      :  Bio.Reactamole.Examples
-- Copyright   :  (c) TBD
-- License     :  TBD
--
-- Maintainer  :  TBD
-- Stability   :  TBD
-- Portability :  TBD
--
-- Examples of uses for Reactamole.
--
--  * Including all examples presented in:
--    /Reactamole: Functional Reactive Molecular Programming/, by Kinge,
--    Lathrop, Osera and Rogers in /DNA Computing and Molecular Programming/ 27,
--    September 2021.

module Bio.Reactamole.Examples
  ( -- * Example SFs
    sinSF
  , lowPass
  , bandPass
  , modulate
  , demodulate
  , srLatch
  , clock
  , unanimousSF

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

-- | SF for the sine function.
sinSF :: SF a Double
sinSF = loop $ sndSF >>> negateSF >>> integrate 1 >>> integrate 0 >>> dupSF

-- | Multiply a 'Signal' by the given constant.
constMult :: Double -> SF Double Double
constMult d = constRl d &&& idSF >>> multSF

-- lowPass1 :: SF Double Double
-- lowPass1 = loop (second negateSF >>> addSF >>> intSF >>> dupSF)

-- | A lowpass filter.
lowPass :: Double -> Double -> SF Double Double
lowPass a b = loop (constMult a *** constMult (-b) >>> addSF >>> integrate 0 >>> dupSF)


-- lowPass2 :: Double -> Double -> SF Double Double
-- lowPass2 k w = loop (constMult k *** constMult (-w) >>> addSF >>> intSF >>> dupSF)

-- bandPass1 :: SF Double Double
-- bandPass1 = loop (second (negateSF &&& intSF >>> addSF) >>> addSF >>> intSF >>> dupSF)

-- bandPass2 :: Double -> Double -> Double -> SF Double Double
-- bandPass2 k q w =
--   loop (second (negateSF &&& intSF >>> (constMult w *** constMult q >>> addSF))
--         >>> first (constMult k) >>> addSF >>> intSF >>> dupSF)

-- bandPass3 :: Double -> Double -> Double -> SF Double Double
-- bandPass3 k q w = loop (first (constMult k)
--   >>> second (constMult (-w) &&& (intSF >>> constMult q) >>> addSF)
--   >>> addSF >>> intSF >>> dupSF)

-- | A bandpass filter.
bandPass :: Double -> Double -> Double -> SF Double Double
bandPass a b c = loop (first (constMult a)
  >>> second (constMult (-c) &&& (integrate 0 >>> constMult (-b)) >>> addSF)
  >>> addSF >>> integrate 0 >>> dupSF)

-- | Generate a sine wave with the given frequency.
carrier :: Double -> SF a Double
carrier w = loop (sndSF >>> constMult w >>> integrate 1
        >>> constMult (-w) >>> integrate 0
        >>> dupSF)

-- | Produce a pair @m(t)@ that satisfies
-- \( \frac{dm}{dt} = u(t)\cdot\sin(ft) - m(t) \).
modulate :: Double -> SF Double Double
modulate w = loop (first (idSF &&& carrier w >>> multSF)
  >>> second negateSF >>> addSF >>> integrate 0 >>> dupSF)

-- | Demodulate a signal that has been modulated on a carrier signal at
-- frequency @w@.
demodulate :: Double -> Double -> SF Double Double
demodulate w q = bandPass (w/q) (w/q) (w*w) >>> rectify >>> lowPass w w

-- expSF :: SF a Double
-- expSF = loop (sndSF >>> intSF >>> dupSF)

-- tanh :: SF a Double
-- tanh = loop (sndSF >>> square >>> negateSF >>> plusOne >>> intSF >>> dupSF)
--   where square = dupSF >>> multSF
--         plusOne = constRl 1 &&& idSF >>> addSF

-- bandPass4 :: Double -> Double -> Double -> SF Double Double
-- bandPass4 a b c =
--   loop
--     ( first (constMult a)
--         >>> second (constMult (- c) &&& (intSF >>> constMult (- b)) >>> addSF)
--         >>> addSF
--         >>> intSF
--         >>> dupSF
--     )

rectify :: SF Double Double
rectify = isPos &&& dupSF >>> entangle >>> (idSF ||| constRl 0)

-- | Determine if input double is positive.
isPos :: SF Double Bool
isPos = posSF

-- boolClock :: Double -> SF () Bool
-- boolClock t = sin  >>> isPos

--------------------------------------------------------------------------------

-- | An SR latch SF.
srLatch :: SF (Bool, Bool) (Bool, Bool)
srLatch = loop $
  arrSp (\(PairS (PairS s' r') (PairS q q'))
         -> PairS (PairS s' q') (PairS r' q))
  >>> (nandSF *** nandSF)
  >>> dupSF

-- | A clock SF.
clock :: SF a Bool
clock = sinSF >>> isPos

-- | A Haskell function to determine if three inputs are the same.
--
-- For demonstration of lifting Boolean functions.
unanimous :: Bool -> Bool -> Bool -> Bool
unanimous x y z = if x then y && z else not (y || z)

-- | An SF to determine if three inputs are the same.
--
-- Demonstrates the use of arr3Bl.
unanimousSF :: SF (Bool, Bool, Bool) Bool
unanimousSF = arr3Bl unanimous
