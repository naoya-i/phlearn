
; Example from [Jensen 96]
(B (name ax1) (=> (^ (rainy x) (w x) ) (watson_wet x)) )
(B (name ax2) (=> (sprinkler x :1.2) (holmes_wet x)) )
(B (name ax3) (=> (sprinkler x :1.2) (kids_play x)) )
;(B (name ax3) (=> (rainy x :1.2) (watson_wet x)) )
(B (name ax4) (=> (rainy x :0.8) (holmes_wet x)) )

(O (name tester01) (^ (holmes_wet A) (kids_play A) (watson_wet A) ) )
;(O (name tester01) (^ (r A :10) ) )

