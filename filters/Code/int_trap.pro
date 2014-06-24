FUNCTION int_trap, x, y

; IDL> x=[1,2,4,6,8]
; IDL> y=x^2

; trapizium rule integration for the case
; where x intervals may not the same

on_error,2



if n_elements(x) eq 0 then begin
    int_trap=0
    message,'scalar input',/inf
endif

; sorting arrays

order = sort(x)

x0 = x(order)
y0 = y(order)


; differencing

dx = -ts_diff(x0,1)
;dy = -ts_diff(y0,1)

; 

; y value at the next x point
yx_1 = shift(y0, -1)

; delta area = dx0*y0 +  1/2 *dx0*dy = dx0/2(y0+y1)

; N.B. we can do the total of all elements as the last one
; is set to zero since dx [ nmax] =0

int_trap = total(dx/2.*(y0+yx_1))

return, int_trap


END

