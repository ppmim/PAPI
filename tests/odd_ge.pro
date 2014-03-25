function odd_ge,x
   if floor(x) MOD 2 eq 1 then return,floor(x)
   return, odd_ge(floor(x)+1)
end
