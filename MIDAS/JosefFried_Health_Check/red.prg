!
! red.prg 
!


! July 6
goto M51

@@ prepflat jul 59,61 68,70 flat_b6 dark6
@@ prepflat jul 56,58 68,70 flat_v6 dark6
@@ prepflat jul 62,64 68,70 flat_r6 dark6
@@ prepflat jul 65,67 68,70 flat_i6 dark6



@@ fcor jul 75,76 flat_v6 dark6 M51v
@@ fcor jul 77,77 flat_b6 dark6 M51b
@@ fcor jul 78,79 flat_r6 dark6 M51r
@@ fcor jul 80,80 flat_i6 dark6 M51i

@@ fcor jul 81,84 flat_v6 dark6 A2151v
@@ fcor jul 85,88 flat_b6 dark6 A2151b1
@@ fcor jul 93,97 flat_b6 dark6 A2151b2
@@ fcor jul 89,92 flat_i6 dark6 A2151i

@@ fcor jul 98,100 flat_r6 dark6 Pl
@@ fcor jul 34,35 flat_r6 dark6 Pl
@@ fcor jul 39,40 flat_r6 dark6 Pl

@@ fcor jul 101,102 flat_r6 dark6 ngc7332
@@ fcor jul 103,104 flat_r6 dark6 ngc7457
@@ fcor jul 105,106 flat_r6 dark6 ngc7600




comp M51v = M51v0075_1+M51v0076_1
comp M51b = M51b0077_1
comp M51r = M51r0078_1+M51r0079_1
comp M51i = M51b0080_1
comp M51il = log10(M51i)
comp M51rl = log10(M51r)
comp M51vl = log10(M51v)





fil/smoo M51v M51vs 8,8
fil/smoo M51b M51bs 8,8
fil/smoo M51r M51rs 8,8
fil/smoo M51i M51is 8,8
comp M51il = log10(M51is)
comp M51rl = log10(M51rs)
comp M51vl = log10(M51vs)
comp M51bl = log10(M51bs)


@@ ali2 M51rl M51vl M51rlr 
@@ ali2 M51il M51vl M51ilr






fil/smoo M51rlr M51rlrs 150,150
fil/smoo M51vl M51vls 150,150
fil/smoo M51ilr M51ilrs 150,150

comp M51v1 = M51vl-M51vls
comp M51r1 = M51rlr-M51rlrs
comp M51i1 = M51ilr-M51ilrs



fil/med M51vl M51vlm 24,24
fil/med M51rlr M51rlm 24,24
fil/med M51ilr M51ilm 24,24
com M51v2 = M51vl-M51vlm
com M51r2 = M51rlr-M51rlm
com M51i2 = M51ilr-M51ilm




M51:
@@ rgb M51ilr M51rlr M51vl M51a 3,4.2 3.3,4.2 3.768,4.2
@@ rgb M51ilr M51rlr M51vl M51b 3,4   3.3,4.2 3.768,4.2



exit:
