
; Example from [Jensen 96]
(B (name ax1) (=> (^ (rainy x) (w x) ) (watson_wet x)) )
(B (name ax2) (=> (sprinkler x :1.2) (holmes_wet x)) )
(B (name ax3) (=> (sprinkler x :1.2) (kids_play x)) )
;(B (name ax3) (=> (rainy x :1.2) (watson_wet x)) )
(B (name ax4) (=> (rainy x :0.8) (holmes_wet x)) )
;(B (name ax5) (=> (^ (god_cry x :0.8) (foo x)) (rainy x)) )

(O (name tester01) (^ (holmes_wet h) (kids_play k) (watson_wet w) ) )
; (O (name tester01) (^ (r A :10) ) )

(B (xor (c1 x) (c2 x) ))

(B (name clsax1) (=> (c1 x) (f1 x)))
(B (name clsax2) (=> (c1 x) (f2 x)))
(B (name clsax3) (=> (c1 x) (f3 x)))
(B (name clsax4) (=> (c2 x) (f1 x)))
(B (name clsax5) (=> (c2 x) (f2 x)))
(B (name clsax6) (=> (c2 x) (f3 x)))

(O (name cls1) (^ (f1 A) (f2 A) (f2 B) (f3 B) ) )
(O (name cls2) (^ (f1 A) (f2 A) (f1 B) (f2 B) ) )
(O (name cls3) (^ (f2 A) (f3 A) (f1 B) (f2 B) ) )
(O (name cls4) (^ (f3 A) (f2 A) ) )

(O (name clsm1) (^ (f1 A) (f2 A) (f2 B) (f3 B) ) )
(O (name clsm2) (^ (f1 A) (f2 A) (f1 B) (f2 B) ) )
(O (name clsm3) (^ (f2 A) (f3 A) (f1 B) (f2 B) ) )
(O (name clsm4) (^ (f3 A) (f2 A) ) )

