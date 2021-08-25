FUNCTION gaussians_sum, X, A
 RETURN, [[A[0]*exp(-((X-A[1])/A[2])^2.0)], [exp(-((X-A[1])/A[2])^2.0)], $
       [2.0*A[0]*exp(-((X-A[1])/A[2])^2.0)*((X-A[1])/(A[2])^2.0)], 
       [2.0*A[0]*exp(-((X-A[1])/A[2])^2.0)*((X-A[1])^2.0/(A[2])^3.0)]]
END

; A = [Amplitude, mi, sigma]
; derivada parcial da gaussiano em relacao A[0]
; derivada parcial da gaussiano em relacao A[1]
; derivada parcial da gaussiano em relacao A[2]

; FUNCTION file_lines, filename
;    OPENR, unit, filename, /GET_LUN
;    str = ''
;    count = 0ll
;    WHILE ~ EOF(unit) DO BEGIN
;       READF, unit, str
;       count = count + 1
;    ENDWHILE
;    FREE_LUN, unit
;    RETURN, count
; END

; Please cite Schlawin et al. (2010) if you make use of this IDL code
;

; TABLE II, pag. 5 -> Chebyshev Approximations for the Complete Elliptic Integrals K and E 
; by W. J. Cody

function ellk,k
; Computes polynomial approximation for the complete elliptic
; integral of the first kind (Hasting's approximation):
m1=1.d0-k^2
a0=1.38629436112d0   ; = ln(4)
a1=0.09666344259d07
a2=0.03590092383d0
a3=0.03742563713d0
a4=0.01451196212d0
b0=0.5d0
b1=0.12498593597d0
b2=0.06880248576d0
b3=0.03328355346d0
b4=0.00441787012d0
ek1=a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
ek2=(b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*alog(m1)
return,ek1-ek2
end

function ellec,k
; Computes polynomial approximation for the complete elliptic
; integral of the second kind (Hasting's approximation):
m1=1.d0-k^2
a1=0.44325141463d0
a2=0.06260601220d0
a3=0.04757383546d0
a4=0.01736506451d0
b1=0.24998368310d0
b2=0.09200180037d0
b3=0.04069697526d0
b4=0.00526449639d0
ee1=1.d0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
ee2=m1*(b1+m1*(b2+m1*(b3+m1*b4)))*alog(1.d0/m1)
return,ee1+ee2
end

function ellpic_bulirsch,n,k
;+
; NAME:
;   ELLPIC_BULIRSCH
;
; PURPOSE:
;   Computes the complete elliptical integral of the third kind using
;   the algorithm of Bulirsch (1965):
;
;   Bulirsch 1965, Numerische Mathematik, 7, 78
;   Bulirsch 1965, Numerische Mathematik, 7, 353
;
; CALLING SEQUENCE:
;    result = ellpic_bulirsch(n, k)
;
; INPUTS:
;
;    n,k - int(dtheta/((1-n*sin(theta)^2)*sqrt(1-k^2*sin(theta)^2)),0, pi/2)
;
; RESULT:
;
;    The complete elliptical integral of the third kind
;
; MODIFICATION HISTORY
;
;  2009/03 -- Written by Eric Agol
;-

kc=sqrt(1d0-k^2) &  p=n+1d0
if(min(p) lt 0d0) then print,'Negative p'
m0=1d0 & c=1d0 & p=sqrt(p) & d=1d0/p & e=kc
niter=0
while niter lt 20 do begin
    f = c & c = d/p+c & g = e/p & d = 2d0*(f*g+d)
    p = g + p & g = m0 & m0 = kc + m0
    if(max(abs(1d0-kc/g)) gt 1.d-8) then begin
        kc = 2d0*sqrt(e) & e=kc*m0
    endif else return,0.5d0*!dpi*(c*m0+d)/(m0*(m0+p))
    niter++
endwhile
return,0.5d0*!dpi*(c*m0+d)/(m0*(m0+p))
end

pro simula_curvas_luz_plus

;nome (caminho completo) da curva de luz correspondente ao eclipse m�dio, autocalibrado
curva_observada='C:\RSI\IDL61\dados_Corot\Curvas_luz_exoplanetas\ID100725706\Butterworth_n2_f02\curva_luz_eclipse_medio_ID100725706_Butterworth_n2_f02_autocalibrada.txt'

;intervalos de amostragem dos valores dos par�metros
delta_b=0.01
delta_adivR=0.1
delta_periodo=0.01
delta_p=0.01

;nomes (caminhos completos) dos arquivos de par�metros

;valores de x (coordenada do planeta, ao longo do eixo x, em fun��o do raio da estrela)
valores_x='C:\RSI\IDL61\dados_Corot\Curvas_luz_exoplanetas\ID100725706\Butterworth_n2_f02\Valores_x_simulacao.txt'

;valores do par�metro de impacto do tr�nsito
valores_b='C:\RSI\IDL61\dados_Corot\Curvas_luz_exoplanetas\ID100725706\Butterworth_n2_f02\Valores_b_simulacao_rodada2.txt'

;valores do raio do planeta em compara��o com o da estrela
valores_p='C:\RSI\IDL61\dados_Corot\Curvas_luz_exoplanetas\ID100725706\Butterworth_n2_f02\Valores_p_simulacao_rodada2.txt'

;valores do per�odo orbital a serem considerados
valores_periodo='C:\RSI\IDL61\dados_Corot\Curvas_luz_exoplanetas\ID100725706\Butterworth_n2_f02\Valores_periodo_simulacao_rodada2.txt'

;valores do raio orbital em compara��o com o raio da estrela
valores_adivR='C:\RSI\IDL61\dados_Corot\Curvas_luz_exoplanetas\ID100725706\Butterworth_n2_f02\Valores_adivR_simulacao_rodada2.txt'

;nome (caminho completo) da curva simulada a ser gerada correspondente ao melhor modelo
curva_saida='C:\RSI\IDL61\dados_Corot\Curvas_luz_exoplanetas\ID100725706\Butterworth_n2_f02\curva_simulada_ID100725706_medio_Butterworth_n2_f02_rodada3_Guilherme.txt'

;nomes (caminhos completos) das tabelas a serem geradas
tabela_saida='C:\RSI\IDL61\dados_Corot\Curvas_luz_exoplanetas\ID100725706\Butterworth_n2_f02\tabela_quis_quadrados_ID100725706_medio_Butterworth_n2_f02_rodada3_Guilherme.txt'
tabela_saida_ordenada='C:\RSI\IDL61\dados_Corot\Curvas_luz_exoplanetas\ID100725706\Butterworth_n2_f02\tabela_quis_quadrados_ID100725706_medio_Butterworth_n2_f02_ordenada_rodada3_Guilherme.txt'
tabela_saida_resultados='C:\RSI\IDL61\dados_Corot\Curvas_luz_exoplanetas\ID100725706\Butterworth_n2_f02\tabela_resultados_ID100725706_medio_Butterworth_n2_f02_rodada3_Guilherme.txt'

;nomes (caminhos completos) dos histogramas a serem gerados
saida_histograma_adivR='C:\RSI\IDL61\dados_Corot\Curvas_luz_exoplanetas\ID100725706\Butterworth_n2_f02\histograma_adivR_ID100725706_medio_Butterworth_n2_f02_rodada3_Guilherme.txt'
saida_histograma_b='C:\RSI\IDL61\dados_Corot\Curvas_luz_exoplanetas\ID100725706\Butterworth_n2_f02\histograma_b_ID100725706_medio_Butterworth_n2_f02_rodada3_Guilherme.txt'
saida_histograma_p='C:\RSI\IDL61\dados_Corot\Curvas_luz_exoplanetas\ID100725706\Butterworth_n2_f02\histograma_p_ID100725706_medio_Butterworth_n2_f02_rodada3_Guilherme.txt'
saida_histograma_periodo='C:\RSI\IDL61\dados_Corot\Curvas_luz_exoplanetas\ID100725706\Butterworth_n2_f02\histograma_periodo_ID100725706_medio_Butterworth_n2_f02_rodada3_Guilherme.txt'

;coeficientes do limb darkening (a princ�pio, n�o mexer aqui)
gamma1=0.44
gamma2=0.23


;in�cio do programa
Nx=file_lines(valores_x)
xs=make_array(Nx,/FLOAT)
Nb=file_lines(valores_b)
bs=make_array(Nb,/FLOAT)
Np=file_lines(valores_p)
ps=make_array(Np,/FLOAT)
Nperiodo=file_lines(valores_periodo)
periodos=make_array(Nperiodo,/FLOAT)
NadivR=file_lines(valores_adivR)
adivRs=make_array(NadivR,/FLOAT)
Nobs=file_lines(curva_observada)
curva_observada_eclipse=make_array(3,Nobs,/double)
tabela_final=make_array(5,/FLOAT)


curva_simulada=make_array(2,Nx,/FLOAT)
curva_simulada_reamostrada=make_array(2,Nobs,/FLOAT)

F=make_array(Nx,/FLOAT)
z=make_array(Nx,/FLOAT)

OPENR,1,curva_observada
READF,1,curva_observada_eclipse
close,1

OPENR,1,valores_x
READF,1,xs
close,1

OPENR,1,valores_b
READF,1,bs
close,1

OPENR,1,valores_p
READF,1,ps
close,1

OPENR,1,valores_periodo
READF,1,periodos
close,1

OPENR,1,valores_adivR
READF,1,adivRs
close,1

OPENW,1,tabela_saida


;simulação das curvas (com limb darkening) e cálculo dos quis quadrados
v=0.0
FOR i=0, Nb-1 DO BEGIN
   b_impacto=bs[i]
   ;print, b
   z=SQRT(xs^2.0+b_impacto^2.0)
   FOR j=0, Np-1 DO BEGIN
      p=ps[j]
      FOR l=0, Nperiodo-1 DO BEGIN
         periodo=periodos[l]
         FOR m=0, NadivR-1 DO BEGIN
            adivR=adivRs[m]
            FOR w=0, Nx-1 DO BEGIN
               ;aplicaçao da modelagem básica
               IF (1.0+p LT z[w]) THEN BEGIN
                  lambda_e=0.0
               ENDIF
               IF (ABS(1.0-p) LT z[w]) and (z[w] LE 1.0+p) THEN BEGIN
                  k0=ACOS((p^2.0+z[w]^2.0-1.0)/(2.0*p*z[w]))
                  k1=ACOS((1.0-p^2.0+z[w]^2.0)/(2.0*z[w]))
                  lambda_e=(1.0/!PI)*(p^2.0*k0+k1-SQRT((4.0*z[w]^2.0-(1.0+z[w]^2.0-p^2.0)^2.0)/4.0))
               ENDIF
               IF (z[w] LE 1.0-p) THEN BEGIN
                  lambda_e=p^2.0
               ENDIF
               IF (z[w] LE p-1.0) THEN BEGIN
                  lambda_e=1.0
               ENDIF


               ;aplicação do limb darkening
               a=(z[w]-p)^2.0
               b=(z[w]+p)^2.0
               k=SQRT((1.0-a)/(4.0*z[w]*p))
               q=p^2.0-z[w]^2.0

               c1=0.0
               c3=0.0
               c2=gamma1+2.0*gamma2
               c4=-1.0*gamma2
               c0=1.0-c1-c2-c3-c4
               Omega=(c0/(4.0))+(c1/(5.0))+(c2/(6.0))+(c3/(7.0))+(c4/(8.0))
               IF (p GT 0.0) and (z[w] GE 1.0+p) THEN BEGIN
                  lambda_d=0.0
                  eta_d=0.0
                  IF (p LE z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
                  ENDIF
                  IF (p GT z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
                  ENDIF
               ENDIF
               IF (p EQ 0.0) and (z[w] GE 0.0) THEN BEGIN
                  lambda_d=0.0
                  eta_d=0.0
                  IF (p LE z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
                  ENDIF
                  IF (p GT z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
                  ENDIF
               ENDIF

               IF (p GT 0.0) and (z[w] GT 0.5+ABS(p-0.5)) and (z[w] LT 1.0+p) THEN BEGIN
                  ;print,(a-1.0)/a
                  ;print, k
                  ;print, ellpic_bulirsch(((a-1.0)/a),k)
                  lambda1=(1.0/(9.0*!PI*SQRT(p*z[w]))) * (((1.0-b)*(2.0*b+a-3.0)-3.0*q*(b-2.0))*ellk(k)+4.0*p*z[w]*(z[w]^2.0+7.0*p^2.0-4.0)*ellec(k)-(3.0*q/a)*ellpic_bulirsch(ABS((a-1.0)/a),k))
                  eta_2=(p^2.0/2.0)*(p^2.0+2.0*z[w]^2.0)
                  eta_1=(2.0*!PI)^(-1.0)*(k1+2.0*eta_2*k0-0.25*(1.0+5.0*p^2.0+z[w]^2.0)*SQRT((1.0-a)*(b-1.0)))
                  lambda_d=lambda1
                  eta_d=eta_1
                  IF (p LE z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
                  ENDIF
                  IF (p GT z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
                  ENDIF
               ENDIF

               ;case III
               IF (p GT 0.0) and (p LT 0.5) and (z[w] GT p) and (z[w] LT 1.0-p) THEN BEGIN
                  lambda2=(2.0/(9.0*!PI*SQRT(1.0-a)))*((1.0-5.0*z[w]^2.0+p^2.0+q^2.0)*ellk(k^(-1.0))+(1.0-a)*(z[w]^2.0+7.0*p^2.0-4.0)*ellec(k^(-1.0))-(3.0*q/a)* ellpic_bulirsch(ABS((a-b)/a),k^(-1.0)))
                  eta_2=(p^2.0/2.0)*(p^2.0+2.0*z[w]^2.0)
                  eta_1=(2.0*!PI)^(-1.0)*(k1+2.0*eta_2*k0-0.25*(1.0+5.0*p^2.0+z[w]^2.0)*SQRT((1.0-a)*(b-1.0)))
                  lambda_d=lambda2
                  eta_d=eta_2
                  IF (p LE z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
                  ENDIF
                  IF (p GT z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
                  ENDIF
               ENDIF

               ;case IV
               IF (p GT 0.0) and (p LT 0.5) and (z[w] EQ 1.0-p) THEN BEGIN
                  ;print, (a-b)/a
                  ;print, k^(-1.0)
                  ;print, ellpic_bulirsch(((a-b)/a),k^(-1.0))
                  IF (p LE 0.5) THEN BEGIN
                     lambda5=(2.0/(3.0*!PI))*ACOS(1.0-2.0*p)-(4.0/(9.0*!PI))*(3.0+2.0*p-8.0*p^2.0)*SQRT(p*(1.0-p))
                  ENDIF
                  IF (p GT 0.5) THEN BEGIN
                     lambda5=(2.0/(3.0*!PI))*ACOS(1.0-2.0*p)-(4.0/(9.0*!PI))*(3.0+2.0*p-8.0*p^2.0)*SQRT(p*(1.0-p))-(2.0/3.0)
                  ENDIF
                  eta_2=(p^2.0/2.0)*(p^2.0+2.0*z[w]^2.0)
                  eta_1=(2.0*!PI)^(-1.0)*(k1+2.0*eta_2*k0-0.25*(1.0+5.0*p^2.0+z[w]^2.0)*SQRT((1.0-a)*(b-1.0)))
                  lambda_d=lambda5
                  eta_d=eta_2
                  IF (p LE z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
                  ENDIF
                  IF (p GT z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
                  ENDIF
               ENDIF

               ;case V
               IF (p GT 0.0) and (p LT 0.5) and (z[w] EQ p) THEN BEGIN
                  ;print, k^(-1.0)
                  ;print, ellpic_bulirsch(((a-b)/a),k^(-1.0))
                  lambda4=(1.0/3.0)+(2.0/(9.0*!PI))*(4.0*(2.0*p^2.0-1.0)*ellec(2.0*p)+(1.0-4.0*p^2.0)*ellk(2.0*p))
                  eta_2=(p^2.0/2.0)*(p^2.0+2.0*z[w]^2.0)
                  eta_1=(2.0*!PI)^(-1.0)*(k1+2.0*eta_2*k0-0.25*(1.0+5.0*p^2.0+z[w]^2.0)*SQRT((1.0-a)*(b-1.0)))
                  lambda_d=lambda4
                  eta_d=eta_2
                  IF (p LE z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
                  ENDIF
                  IF (p GT z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
                  ENDIF
               ENDIF

               ;case VI
               IF (p EQ 0.5) and (z[w] EQ 0.5) THEN BEGIN
                  ;print, (a-b)/a
                  ;print, k^(-1.0)
                  ;print, ellpic_bulirsch(((a-b)/a),k^(-1.0))
                  lambda_d=(1.0/3.0)-(4.0/(9.0*!PI))
                  eta_d=3.0/32.0
                  IF (p LE z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
                  ENDIF
                  IF (p GT z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
                  ENDIF
               ENDIF

               ;case VII
               IF (p GT 0.5) and (z[w] EQ p) THEN BEGIN
                  ;print, (a-b)/a
                  ;print, k^(-1.0)
                  ;print, ellpic_bulirsch(((a-b)/a),k^(-1.0))
                  lambda3=(1.0/3.0)+((16.0*p)/(9.0*!PI))*(2.0*p^2.0-1.0)*ellec(1.0/(2.0*p))-(((1.0-4.0*p^2.0)*(3.0-8.0*p^2.0))/(9.0*!PI*p))*ellk(1.0/(2.0*p))
                  eta_2=(p^2.0/2.0)*(p^2.0+2.0*z[w]^2.0)
                  eta_1=(2.0*!PI)^(-1.0)*(k1+2.0*eta_2*k0-0.25*(1.0+5.0*p^2.0+z[w]^2.0)*SQRT((1.0-a)*(b-1.0)))
                  lambda_d=lambda3
                  eta_d=eta_1
                  IF (p LE z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
                  ENDIF
                  IF (p GT z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
                  ENDIF
               ENDIF

               ;case VIII
               IF (p GT 0.5) and (z[w] GE ABS(1.0-p)) and (z[w] LT p) THEN BEGIN
                  ;print, (a-b)/a
                  ;print, k^(-1.0)
                  ;print, ellpic_bulirsch(((a-b)/a),k^(-1.0))
                  lambda1=(1.0/(9.0*!PI*SQRT(p*z[w])))*(((1.0-b)*(2.0*b+a-3.0)-3.0*q*(b-2.0))*ellk(k)+4.0*p*z[w]*(z[w]^2.0+7.0*p^2.0-4.0)*ellec(k)-(3.0*q/a)*ellpic_bulirsch(ABS((a-1.0)/a),k))
                  eta_2=(p^2.0/2.0)*(p^2.0+2.0*z[w]^2.0)
                  eta_1=(2.0*!PI)^(-1.0)*(k1+2.0*eta_2*k0-0.25*(1.0+5.0*p^2.0+z[w]^2.0)*SQRT((1.0-a)*(b-1.0)))
                  lambda_d=lambda1
                  eta_d=eta_1
                  IF (p LE z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
                  ENDIF
                  IF (p GT z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
                  ENDIF
               ENDIF

               ;case IX
               IF (p GT 0.0) and (p LT 1.0) and (z[w] GT 0.0) and (z[w] LT 0.5-ABS(p-0.5)) THEN BEGIN
                  ;print, (a-b)/a
                  ;print, k^(-1.0)
                  ;print, ellpic_bulirsch(((a-b)/a),k^(-1.0))
                  lambda2=(2.0/(9.0*!PI*SQRT(1.0-a)))*((1.0-5.0*z[w]^2.0+p^2.0+q^2.0)*ellk(k^(-1.0))+(1.0-a)*(z[w]^2.0+7.0*p^2.0-4.0)*ellec(k^(-1.0))-(3.0*q/a)* ellpic_bulirsch(ABS((a-b)/a),k^(-1.0)))
                  eta_2=(p^2.0/2.0)*(p^2.0+2.0*z[w]^2.0)
                  eta_1=(2.0*!PI)^(-1.0)*(k1+2.0*eta_2*k0-0.25*(1.0+5.0*p^2.0+z[w]^2.0)*SQRT((1.0-a)*(b-1.0)))
                  lambda_d=lambda2
                  eta_d=eta_2
                  IF (p LE z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
                  ENDIF
                  IF (p GT z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
                  ENDIF
               ENDIF

               ;case X
               IF (p GT 0.0) and (p LT 1.0) and (z[w] EQ 0.0) THEN BEGIN
                  ;print, (a-b)/a
                  ;print, k^(-1.0)
                  ;print, ellpic_bulirsch(((a-b)/a),k^(-1.0))
                  lambda6=-(2.0/3.0)*(1.0-p^2.0)^(3.0/2.0)
                  eta_2=(p^2.0/2.0)*(p^2.0+2.0*z[w]^2.0)
                  eta_1=(2.0*!PI)^(-1.0)*(k1+2.0*eta_2*k0-0.25*(1.0+5.0*p^2.0+z[w]^2.0)*SQRT((1.0-a)*(b-1.0)))
                  lambda_d=lambda6
                  eta_d=eta_2
                  IF (p LE z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
                  ENDIF
                  IF (p GT z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
                  ENDIF
               ENDIF

               ;case XI
               IF (p GT 1.0) and (z[w] GE 0.0) and (z[w] LT p-1.0) THEN BEGIN
                  ;print, (a-b)/a
                  ;print, k^(-1.0)
                  ;print, ellpic_bulirsch(((a-b)/a),k^(-1.0))
                  lambda_d=0.0
                  eta_d=0.5
                  IF (p LE z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
                  ENDIF
                  IF (p GT z[w]) THEN BEGIN
                     F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
                  ENDIF
               ENDIF
            ENDFOR
            ;;;

            ;conversao de x para tempo
            ttrans=periodo/(!PI*adivR)
            meio = (Nx/2.0)-0.5
            FOR u=0, Nx-1 DO BEGIN
               IF (u LT meio) THEN BEGIN
                  curva_simulada[0,u]=(-1.0)*xs[u]*ttrans/2.0
               ENDIF
               IF (u GE meio) THEN BEGIN
                  curva_simulada[0,u]=xs[u]*ttrans/2.0
               ENDIF
            ENDFOR
            curva_simulada[1,*]=F[*]

            ;reamostragem
            curva_simulada_reamostrada[0,*]=curva_observada_eclipse[0,*]
            curva_simulada_reamostrada[1,*]=INTERPOL(curva_simulada[1,*],curva_simulada[0,*],curva_simulada_reamostrada[0,*])


            ;calculo do qui2
            qui2=0.0
            FOR u=0, Nobs-1 DO BEGIN
               qui2=qui2+((curva_observada_eclipse[1,u]-curva_simulada_reamostrada[1,u])^2.0/curva_observada_eclipse[2,u]^2.0)
            ENDFOR

            ;print, curva_observada_eclipse[0,*]
            tabela_final[0]=b_impacto
            tabela_final[1]=p
            tabela_final[2]=periodo
            tabela_final[3]=adivR
            tabela_final[4]=qui2
            string_v=STRING(v+1.0)
            string_v=STRTRIM(string_v,1)
            v=v+1.0
            PRINTF,1,tabela_final[*]
         ENDFOR
      ENDFOR
   ENDFOR
ENDFOR

close,1

;ordenamento da tabela de quis quadrados
Nmax=float(Nb)*float(Nperiodo)*float(Np)*float(NadivR)
nova_tabela_final=make_array(5,Nmax,/FLOAT,value=0.0)
nova_tabela_final_ordenada=make_array(5,Nmax,/FLOAT,value=0.0)

OPENR,1,tabela_saida
READF,1,nova_tabela_final
close,1

FOR n=0.0, Nmax-1 DO BEGIN
   result=MAX(nova_tabela_final[4,*],Max_subscript)
   nova_tabela_final_ordenada[*,Nmax-1-n]=nova_tabela_final[*,Max_subscript]
   nova_tabela_final[*,Max_subscript]=0.0
ENDFOR

n=0.0
WHILE (nova_tabela_final_ordenada[4,n] EQ 0.0) DO BEGIN
   n=n+1.0
ENDWHILE

new_Nmax=(Nmax-1.0)-n+1
nova_tabela_final_ordenada_mesmo=make_array(5,new_Nmax,/FLOAT,value=0.0)

nova_tabela_final_ordenada_mesmo[*,*]=nova_tabela_final_ordenada[*,n:Nmax-1]

n=0.0
WHILE (nova_tabela_final_ordenada_mesmo[4,n]-nova_tabela_final_ordenada_mesmo[4,0] LE 1.0) DO BEGIN
   n=n+1.0
ENDWHILE

newnew_Nmax=n-1.0

OPENW,1,tabela_saida_ordenada
FOR n=0.0, new_Nmax-1 DO BEGIN
   PRINTF,1,nova_tabela_final_ordenada_mesmo[*,n]
ENDFOR
close,1


;confec��o da tabela dos resultados e determina��o das incertezas (atrav�s dos ajustes de histogramas)

;histograma de adivR
tabela_resultados=make_array(2,4,/FLOAT,value=0.0)

adivR_minimo=MIN(nova_tabela_final_ordenada_mesmo[3,0:newnew_Nmax])
adivR_maximo=MAX(nova_tabela_final_ordenada_mesmo[3,0:newnew_Nmax])
numero_adivR=(adivR_maximo-adivR_minimo+delta_adivR)/delta_adivR
IF (fix(numero_adivR) LT numero_adivR) THEN BEGIN
   numero_adivR=fix(numero_adivR)+1
ENDIF
IF (fix(numero_adivR) GT numero_adivR) THEN BEGIN
   numero_adivR=fix(numero_adivR)-1
ENDIF
IF (fix(numero_adivR) EQ numero_adivR) THEN BEGIN
   numero_adivR=fix(numero_adivR)
ENDIF
histograma_final_adivR=make_array(2,numero_adivR,/FLOAT,value=0.0)

FOR n=0.0, numero_adivR-1.0 DO BEGIN
   histograma_final_adivR[0,n]=adivR_minimo+n*delta_adivR
   FOR w=0, newnew_Nmax DO BEGIN
      IF (ABS(nova_tabela_final_ordenada_mesmo[3,w] - histograma_final_adivR[0,n]) LT delta_adivR/10.0) THEN BEGIN
         histograma_final_adivR[1,n]=histograma_final_adivR[1,n]+1.0
      ENDIF
   ENDFOR
ENDFOR

OPENW,1,saida_histograma_adivR
PRINTF,1,histograma_final_adivR
close,1

cont=0.0
FOR n=0.0, numero_adivR-1.0 DO BEGIN
   IF (histograma_final_adivR[1,n] EQ 0.0) THEN BEGIN
      cont=cont+1.0
   ENDIF
ENDFOR

histograma_final_adivR_semzeros=make_array(2,numero_adivR-cont,/FLOAT)

w=0
FOR n=0.0, numero_adivR-1.0 DO BEGIN
   IF (histograma_final_adivR[1,n] NE 0.0) THEN BEGIN
      histograma_final_adivR_semzeros[*,w]=histograma_final_adivR[*,n]
      w=w+1.0
   ENDIF
ENDFOR

print, 'Histograma de adivR, sem zeros:'
print, histograma_final_adivR_semzeros
stop

PLOT,histograma_final_adivR_semzeros[0,*],histograma_final_adivR_semzeros[1,*]

minimo=0.0
maximo=numero_adivR-cont-1
ok=''
READ,ok,PROMPT='Is the range ok?'
WHILE (ok EQ 'no') DO BEGIN
   READ,minimo,PROMPT='Nm�nimo='
   READ,maximo,PROMPT='Nm�ximo='
   PLOT,histograma_final_adivR_semzeros[0,minimo:maximo],histograma_final_adivR_semzeros[1,minimo:maximo]
   READ,ok,PROMPT='Is the range ok?'
ENDWHILE

tabela_resultados[0,3]=nova_tabela_final_ordenada_mesmo[3,0]
IF (maximo-minimo+1 GT 3) THEN BEGIN
   A=make_array(3,/FLOAT)
   FITA=make_array(3,/FLOAT,value=1.0)
   A[0]=MAX(histograma_final_adivR_semzeros[1,minimo:maximo],Max_subscript)
   A[1]=histograma_final_adivR_semzeros[0,Max_subscript]
   A[2]=delta_adivR
   Result_gaussians=LMFIT(histograma_final_adivR_semzeros[0,minimo:maximo], histograma_final_adivR_semzeros[1,minimo:maximo], A, /DOUBLE, FUNCTION_NAME='gaussians_sum', iter=iter, itmin=100, FITA=FITA, measure_errors=measure_errors, sigma=sigma)
   PLOT,histograma_final_adivR_semzeros[0,minimo:maximo],histograma_final_adivR[1,minimo:maximo],PSYM=1
   OPLOT,histograma_final_adivR_semzeros[0,minimo:maximo],Result_gaussians
   ok=''
   READ,ok,PROMPT='Is the fit ok?'
   WHILE (ok EQ 'no') DO BEGIN
      amplitude=0.0
      centro=0.0
      sigma=0.0
      READ,amplitude,PROMPT='A[0]='
      READ,centro,PROMPT='A[1]='
      READ,sigma,PROMPT='A[2]='
      A[0]=amplitude
      A[1]=centro
      A[2]=sigma
      Result_gaussians=LMFIT(histograma_final_adivR_semzeros[0,minimo:maximo], histograma_final_adivR[1,minimo:maximo], A, /DOUBLE, FUNCTION_NAME='gaussians_sum', iter=iter, itmin=100, FITA=FITA, measure_errors=measure_errors, sigma=sigma)
      PLOT,histograma_final_adivR_semzeros[0,minimo:maximo],histograma_final_adivR_semzeros[1,minimo:maximo],PSYM=1
      OPLOT,histograma_final_adivR_semzeros[0,minimo:maximo],Result_gaussians
      READ,ok,PROMPT='Is the fit ok?'
   ENDWHILE
   tabela_resultados[1,3]=A[2]
ENDIF ELSE BEGIN
   tabela_resultados[1,3]=STDDEV(nova_tabela_final_ordenada_mesmo[3,0:newnew_Nmax])
ENDELSE
;;;;



;histograma do per�odo
periodo_minimo=MIN(nova_tabela_final_ordenada_mesmo[2,0:newnew_Nmax])
periodo_maximo=MAX(nova_tabela_final_ordenada_mesmo[2,0:newnew_Nmax])
numero_periodo=(periodo_maximo-periodo_minimo+delta_periodo)/delta_periodo

IF (fix(numero_periodo) LT numero_periodo) THEN BEGIN
   numero_periodo=fix(numero_periodo)+1
ENDIF
IF (fix(numero_periodo) GT numero_periodo) THEN BEGIN
   numero_periodo=fix(numero_periodo)-1
ENDIF
IF (fix(numero_periodo) EQ numero_periodo) THEN BEGIN
   numero_periodo=fix(numero_periodo)
ENDIF
histograma_final_periodo=make_array(2,numero_periodo,/FLOAT,value=0.0)

FOR n=0.0, numero_periodo-1.0 DO BEGIN
   histograma_final_periodo[0,n]=periodo_minimo+n*delta_periodo
   FOR w=0, newnew_Nmax DO BEGIN
      IF (ABS(nova_tabela_final_ordenada_mesmo[2,w] - histograma_final_periodo[0,n]) LT delta_periodo/10.0) THEN BEGIN
         histograma_final_periodo[1,n]=histograma_final_periodo[1,n]+1.0
      ENDIF
   ENDFOR
ENDFOR

OPENW,1,saida_histograma_periodo
PRINTF,1,histograma_final_periodo
close,1


cont=0.0
FOR n=0.0, numero_periodo-1.0 DO BEGIN
   IF (histograma_final_periodo[1,n] EQ 0.0) THEN BEGIN
      cont=cont+1.0
   ENDIF
ENDFOR

histograma_final_periodo_semzeros=make_array(2,numero_periodo-cont,/FLOAT)

w=0
FOR n=0.0, numero_periodo-1.0 DO BEGIN
   IF (histograma_final_periodo[1,n] NE 0.0) THEN BEGIN
      histograma_final_periodo_semzeros[*,w]=histograma_final_periodo[*,n]
      w=w+1.0
   ENDIF
ENDFOR

print, 'Histograma do per�odo, sem zeros:'
print, histograma_final_periodo_semzeros
stop

PLOT,histograma_final_periodo_semzeros[0,*],histograma_final_periodo_semzeros[1,*]


minimo=0.0
maximo=numero_periodo-cont-1
ok=''
READ,ok,PROMPT='Is the range ok?'
WHILE (ok EQ 'no') DO BEGIN
   READ,minimo,PROMPT='Nm�nimo='
   READ,maximo,PROMPT='Nm�ximo='
   PLOT,histograma_final_periodo_semzeros[0,minimo:maximo],histograma_final_periodo_semzeros[1,minimo:maximo]
   READ,ok,PROMPT='Is the range ok?'
ENDWHILE


tabela_resultados[0,2]=nova_tabela_final_ordenada_mesmo[2,0]
IF (maximo-minimo GT 3) THEN BEGIN
   A=make_array(3,/FLOAT)
   FITA=make_array(3,/FLOAT,value=1.0)
   A[0]=MAX(histograma_final_periodo_semzeros[1,minimo:maximo],Max_subscript)
   A[1]=histograma_final_periodo_semzeros[0,Max_subscript]
   A[2]=delta_periodo
   Result_gaussians=LMFIT(histograma_final_periodo_semzeros[0,minimo:maximo], histograma_final_periodo_semzeros[1,minimo:maximo], A, /DOUBLE, FUNCTION_NAME='gaussians_sum', iter=iter, itmin=100, FITA=FITA, measure_errors=measure_errors, sigma=sigma)
   PLOT,histograma_final_periodo_semzeros[0,minimo:maximo],histograma_final_periodo_semzeros[1,minimo:maximo],PSYM=1
   OPLOT,histograma_final_periodo_semzeros[0,minimo:maximo],Result_gaussians
   ok=''
   READ,ok,PROMPT='Is the fit ok?'
   WHILE (ok EQ 'no') DO BEGIN
      amplitude=0.0
      centro=0.0
      sigma=0.0
      READ,amplitude,PROMPT='A[0]='
      READ,centro,PROMPT='A[1]='
      READ,sigma,PROMPT='A[2]='
      A[0]=amplitude
      A[1]=centro
      A[2]=sigma
      Result_gaussians=LMFIT(histograma_final_periodo_semzeros[0,minimo:maximo], histograma_final_periodo_semzeros[1,minimo:maximo], A, /DOUBLE, FUNCTION_NAME='gaussians_sum', iter=iter, itmin=100, FITA=FITA, measure_errors=measure_errors, sigma=sigma)
      PLOT,histograma_final_periodo_semzeros[0,minimo:maximo],histograma_final_periodo_semzeros[1,minimo:maximo],PSYM=1
      OPLOT,histograma_final_periodo_semzeros[0,minimo:maximo],Result_gaussians
      READ,ok,PROMPT='Is the fit ok?'
   ENDWHILE
   tabela_resultados[1,2]=A[2]
ENDIF ELSE BEGIN
   tabela_resultados[1,2]=STDDEV(nova_tabela_final_ordenada_mesmo[2,0:newnew_Nmax])
ENDELSE



;histograma dos valores de p
p_minimo=MIN(nova_tabela_final_ordenada_mesmo[1,0:newnew_Nmax])
p_maximo=MAX(nova_tabela_final_ordenada_mesmo[1,0:newnew_Nmax])
numero_p=(p_maximo-p_minimo+delta_p)/delta_p
IF (fix(numero_p) LT numero_p) THEN BEGIN
   numero_p=fix(numero_p)+1
ENDIF
IF (fix(numero_p) GT numero_p) THEN BEGIN
   numero_p=fix(numero_p)-1
ENDIF
IF (fix(numero_p) EQ numero_p) THEN BEGIN
   numero_p=fix(numero_p)
ENDIF

histograma_final_p=make_array(2,numero_p,/FLOAT,value=0.0)

FOR n=0.0, numero_p-1.0 DO BEGIN
   histograma_final_p[0,n]=p_minimo+n*delta_p
   FOR w=0, newnew_Nmax DO BEGIN
      IF (ABS(nova_tabela_final_ordenada_mesmo[1,w] - histograma_final_p[0,n]) LT delta_p/10.0) THEN BEGIN
         histograma_final_p[1,n]=histograma_final_p[1,n]+1.0
      ENDIF
   ENDFOR
ENDFOR

OPENW,1,saida_histograma_p
PRINTF,1,histograma_final_p
close,1

cont=0.0
FOR n=0.0, numero_p-1.0 DO BEGIN
   IF (histograma_final_p[1,n] EQ 0.0) THEN BEGIN
      cont=cont+1.0
   ENDIF
ENDFOR

histograma_final_p_semzeros=make_array(2,numero_p-cont,/FLOAT)

w=0
FOR n=0.0, numero_p-1.0 DO BEGIN
   IF (histograma_final_p[1,n] NE 0.0) THEN BEGIN
      histograma_final_p_semzeros[*,w]=histograma_final_p[*,n]
      w=w+1.0
   ENDIF
ENDFOR

print, 'Histograma de p, sem zeros:'
print, histograma_final_p_semzeros
stop

PLOT,histograma_final_p_semzeros[0,*],histograma_final_p_semzeros[1,*]


minimo=0.0
maximo=numero_p-cont-1
ok=''
READ,ok,PROMPT='Is the range ok?'
WHILE (ok EQ 'no') DO BEGIN
   READ,minimo,PROMPT='Nm�nimo='
   READ,maximo,PROMPT='Nm�ximo='
   PLOT,histograma_final_p_semzeros[0,minimo:maximo],histograma_final_p_semzeros[1,minimo:maximo]
   READ,ok,PROMPT='Is the range ok?'
ENDWHILE


tabela_resultados[0,1]=nova_tabela_final_ordenada_mesmo[1,0]
IF (maximo-minimo GT 3) THEN BEGIN
   A=make_array(3,/FLOAT)
   FITA=make_array(3,/FLOAT,value=1.0)
   A[0]=MAX(histograma_final_p_semzeros[1,minimo:maximo],Max_subscript)
   A[1]=histograma_final_p_semzeros[0,Max_subscript]
   A[2]=delta_p
   Result_gaussians=LMFIT(histograma_final_p_semzeros[0,minimo:maximo], histograma_final_p_semzeros[1,minimo:maximo], A, /DOUBLE, FUNCTION_NAME='gaussians_sum', iter=iter, itmin=100, FITA=FITA, measure_errors=measure_errors, sigma=sigma)
   PLOT,histograma_final_p_semzeros[0,minimo:maximo],histograma_final_p_semzeros[1,minimo:maximo],PSYM=1
   OPLOT,histograma_final_p_semzeros[0,minimo:maximo],Result_gaussians
   ok=''
   READ,ok,PROMPT='Is the fit ok?'
   WHILE (ok EQ 'no') DO BEGIN
      amplitude=0.0
      centro=0.0
      sigma=0.0
      READ,amplitude,PROMPT='A[0]='
      READ,centro,PROMPT='A[1]='
      READ,sigma,PROMPT='A[2]='
      A[0]=amplitude
      A[1]=centro
      A[2]=sigma
      Result_gaussians=LMFIT(histograma_final_p_semzeros[0,minimo:maximo], histograma_final_p_semzeros[1,minimo:maximo], A, /DOUBLE, FUNCTION_NAME='gaussians_sum', iter=iter, itmin=100, FITA=FITA, measure_errors=measure_errors, sigma=sigma)
      PLOT,histograma_final_p_semzeros[0,minimo:maximo],histograma_final_p_semzeros[1,minimo:maximo],PSYM=1
      OPLOT,histograma_final_p_semzeros[0,minimo:maximo],Result_gaussians
      READ,ok,PROMPT='Is the fit ok?'
   ENDWHILE
   tabela_resultados[1,1]=A[2]
ENDIF ELSE BEGIN
   tabela_resultados[1,1]=STDDEV(nova_tabela_final_ordenada_mesmo[1,0:newnew_Nmax])
ENDELSE
;;;;;;;




;histograma dos valores de b (par�metro de impacto)
b_minimo=MIN(nova_tabela_final_ordenada_mesmo[0,0:newnew_Nmax])
b_maximo=MAX(nova_tabela_final_ordenada_mesmo[0,0:newnew_Nmax])
numero_b=(b_maximo-b_minimo+delta_b)/delta_b

IF (fix(numero_b) LT numero_b) THEN BEGIN
   numero_b=fix(numero_b)+1
ENDIF
IF (fix(numero_b) GT numero_b) THEN BEGIN
   numero_b=fix(numero_b)-1
ENDIF
IF (fix(numero_b) EQ numero_b) THEN BEGIN
   numero_b=fix(numero_b)
ENDIF
histograma_final_b=make_array(2,numero_b,/FLOAT,value=0.0)

FOR n=0.0, numero_b-1.0 DO BEGIN
   histograma_final_b[0,n]=b_minimo+n*delta_b
   FOR w=0, newnew_Nmax DO BEGIN
      IF (ABS(nova_tabela_final_ordenada_mesmo[0,w] - histograma_final_b[0,n]) LT delta_b/10.0) THEN BEGIN
         histograma_final_b[1,n]=histograma_final_b[1,n]+1.0
      ENDIF
   ENDFOR
ENDFOR

OPENW,1,saida_histograma_b
PRINTF,1,histograma_final_b
close,1

cont=0.0
FOR n=0.0, numero_b-1.0 DO BEGIN
   IF (histograma_final_b[1,n] EQ 0.0) THEN BEGIN
      cont=cont+1.0
   ENDIF
ENDFOR

histograma_final_b_semzeros=make_array(2,numero_b-cont,/FLOAT)

w=0
FOR n=0.0, numero_b-1.0 DO BEGIN
   IF (histograma_final_b[1,n] NE 0.0) THEN BEGIN
      histograma_final_b_semzeros[*,w]=histograma_final_b[*,n]
      w=w+1.0
   ENDIF
ENDFOR

print, 'Histograma de b, sem zeros:'
print, histograma_final_b_semzeros
stop

PLOT,histograma_final_b_semzeros[0,*],histograma_final_b_semzeros[1,*]


minimo=0.0
maximo=numero_b-cont-1
ok=''
READ,ok,PROMPT='Is the range ok?'
WHILE (ok EQ 'no') DO BEGIN
   READ,minimo,PROMPT='Nm�nimo='
   READ,maximo,PROMPT='Nm�ximo='
   PLOT,histograma_final_b_semzeros[0,minimo:maximo],histograma_final_b_semzeros[1,minimo:maximo]
   READ,ok,PROMPT='Is the range ok?'
ENDWHILE


tabela_resultados[0,0]=nova_tabela_final_ordenada_mesmo[0,0]
IF (maximo-minimo GT 3) THEN BEGIN
   A=make_array(3,/FLOAT)
   FITA=make_array(3,/FLOAT,value=1.0)
   A[0]=MAX(histograma_final_b_semzeros[1,minimo:maximo],Max_subscript)
   A[1]=histograma_final_b_semzeros[0,Max_subscript]
   A[2]=delta_b
   Result_gaussians=LMFIT(histograma_final_b_semzeros[0,minimo:maximo], histograma_final_b_semzeros[1,minimo:maximo], A, /DOUBLE, FUNCTION_NAME='gaussians_sum', iter=iter, itmin=100, FITA=FITA, measure_errors=measure_errors, sigma=sigma)
   PLOT,histograma_final_b_semzeros[0,minimo:maximo],histograma_final_b_semzeros[1,minimo:maximo],PSYM=1
   OPLOT,histograma_final_b_semzeros[0,minimo:maximo],Result_gaussians
   ok=''
   READ,ok,PROMPT='Is the fit ok?'
   WHILE (ok EQ 'no') DO BEGIN
      amplitude=0.0
      centro=0.0
      sigma=0.0
      READ,amplitude,PROMPT='A[0]='
      READ,centro,PROMPT='A[1]='
      READ,sigma,PROMPT='A[2]='
      A[0]=amplitude
      A[1]=centro
      A[2]=sigma
      Result_gaussians=LMFIT(histograma_final_b_semzeros[0,minimo:maximo], histograma_final_b_semzeros[1,minimo:maximo], A, /DOUBLE, FUNCTION_NAME='gaussians_sum', iter=iter, itmin=100, FITA=FITA, measure_errors=measure_errors, sigma=sigma)
      PLOT,histograma_final_b_semzeros[0,minimo:maximo],histograma_final_b_semzeros[1,minimo:maximo],PSYM=1
      OPLOT,histograma_final_b_semzeros[0,minimo:maximo],Result_gaussians
      READ,ok,PROMPT='Is the fit ok?'
   ENDWHILE
   tabela_resultados[1,0]=A[2]
ENDIF ELSE BEGIN
   tabela_resultados[1,0]=STDDEV(nova_tabela_final_ordenada_mesmo[0,0:newnew_Nmax])
ENDELSE


OPENW,1,tabela_saida_resultados
FOR w=0, 3 DO BEGIN
   PRINTF,1,tabela_resultados[*,w]
ENDFOR
close,1




;elabora��o da curva simulada final
b_impacto=nova_tabela_final_ordenada_mesmo[0,0]
z=SQRT(xs^2.0+b_impacto^2.0)
p=nova_tabela_final_ordenada_mesmo[1,0]
periodo=nova_tabela_final_ordenada_mesmo[2,0]
adivR=nova_tabela_final_ordenada_mesmo[3,0]
FOR w=0, Nx-1 DO BEGIN
   ;aplica�ao da modelagem b�sica
   IF (1.0+p LT z[w]) THEN BEGIN
      lambda_e=0.0
   ENDIF
   IF (ABS(1.0-p) LT z[w]) and (z[w] LE 1.0+p) THEN BEGIN
      k0=ACOS((p^2.0+z[w]^2.0-1.0)/(2.0*p*z[w]))
      k1=ACOS((1.0-p^2.0+z[w]^2.0)/(2.0*z[w]))
      lambda_e=(1.0/!PI)*(p^2.0*k0+k1-SQRT((4.0*z[w]^2.0-(1.0+z[w]^2.0-p^2.0)^2.0)/4.0))
   ENDIF
   IF (z[w] LE 1.0-p) THEN BEGIN
      lambda_e=p^2.0
   ENDIF
   IF (z[w] LE p-1.0) THEN BEGIN
      lambda_e=1.0
   ENDIF


   ;aplica��o do limb darkening
   a=(z[w]-p)^2.0
   b=(z[w]+p)^2.0
   k=SQRT((1.0-a)/(4.0*z[w]*p))
   q=p^2.0-z[w]^2.0

   c1=0.0
   c3=0.0
   c2=gamma1+2.0*gamma2
   c4=-1.0*gamma2
   c0=1.0-c1-c2-c3-c4
   Omega=(c0/(4.0))+(c1/(5.0))+(c2/(6.0))+(c3/(7.0))+(c4/(8.0))
   IF (p GT 0.0) and (z[w] GE 1.0+p) THEN BEGIN
      lambda_d=0.0
      eta_d=0.0
      IF (p LE z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
      ENDIF
      IF (p GT z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
      ENDIF
   ENDIF
   IF (p EQ 0.0) and (z[w] GE 0.0) THEN BEGIN
      lambda_d=0.0
      eta_d=0.0
      IF (p LE z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
      ENDIF
      IF (p GT z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
      ENDIF
   ENDIF

   IF (p GT 0.0) and (z[w] GT 0.5+ABS(p-0.5)) and (z[w] LT 1.0+p) THEN BEGIN
      ;print,(a-1.0)/a
      ;print, k
      ;print, ellpic_bulirsch(((a-1.0)/a),k)
      lambda1=(1.0/(9.0*!PI*SQRT(p*z[w])))*(((1.0-b)*(2.0*b+a-3.0)-3.0*q*(b-2.0))*ellk(k)+4.0*p*z[w]*(z[w]^2.0+7.0*p^2.0-4.0)*ellec(k)-(3.0*q/a)*ellpic_bulirsch(ABS((a-1.0)/a),k))
      eta_2=(p^2.0/2.0)*(p^2.0+2.0*z[w]^2.0)
      eta_1=(2.0*!PI)^(-1.0)*(k1+2.0*eta_2*k0-0.25*(1.0+5.0*p^2.0+z[w]^2.0)*SQRT((1.0-a)*(b-1.0)))
      lambda_d=lambda1
      eta_d=eta_1
      IF (p LE z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
      ENDIF
      IF (p GT z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
      ENDIF
   ENDIF

   ;case III
   IF (p GT 0.0) and (p LT 0.5) and (z[w] GT p) and (z[w] LT 1.0-p) THEN BEGIN
      ;print, (a-b)/a
      ;print, k^(-1.0)
      ;print, ellpic_bulirsch(((a-b)/a),k^(-1.0))
      lambda2=(2.0/(9.0*!PI*SQRT(1.0-a)))*((1.0-5.0*z[w]^2.0+p^2.0+q^2.0)*ellk(k^(-1.0))+(1.0-a)*(z[w]^2.0+7.0*p^2.0-4.0)*ellec(k^(-1.0))-(3.0*q/a)* ellpic_bulirsch(ABS((a-b)/a),k^(-1.0)))
      eta_2=(p^2.0/2.0)*(p^2.0+2.0*z[w]^2.0)
      eta_1=(2.0*!PI)^(-1.0)*(k1+2.0*eta_2*k0-0.25*(1.0+5.0*p^2.0+z[w]^2.0)*SQRT((1.0-a)*(b-1.0)))
      lambda_d=lambda2
      eta_d=eta_2
      IF (p LE z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
      ENDIF
      IF (p GT z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
      ENDIF
   ENDIF

   ;case IV
   IF (p GT 0.0) and (p LT 0.5) and (z[w] EQ 1.0-p) THEN BEGIN
      ;print, (a-b)/a
      ;print, k^(-1.0)
      ;print, ellpic_bulirsch(((a-b)/a),k^(-1.0))
      IF (p LE 0.5) THEN BEGIN
         lambda5=(2.0/(3.0*!PI))*ACOS(1.0-2.0*p)-(4.0/(9.0*!PI))*(3.0+2.0*p-8.0*p^2.0)*SQRT(p*(1.0-p))
      ENDIF
      IF (p GT 0.5) THEN BEGIN
         lambda5=(2.0/(3.0*!PI))*ACOS(1.0-2.0*p)-(4.0/(9.0*!PI))*(3.0+2.0*p-8.0*p^2.0)*SQRT(p*(1.0-p))-(2.0/3.0)
      ENDIF
      eta_2=(p^2.0/2.0)*(p^2.0+2.0*z[w]^2.0)
      eta_1=(2.0*!PI)^(-1.0)*(k1+2.0*eta_2*k0-0.25*(1.0+5.0*p^2.0+z[w]^2.0)*SQRT((1.0-a)*(b-1.0)))
      lambda_d=lambda5
      eta_d=eta_2
      IF (p LE z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
      ENDIF
      IF (p GT z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
      ENDIF
   ENDIF

   ;case V
   IF (p GT 0.0) and (p LT 0.5) and (z[w] EQ p) THEN BEGIN
      ;print, 'aqui'
      ;print, k^(-1.0)
      ;print, ellpic_bulirsch(((a-b)/a),k^(-1.0))
      lambda4=(1.0/3.0)+(2.0/(9.0*!PI))*(4.0*(2.0*p^2.0-1.0)*ellec(2.0*p)+(1.0-4.0*p^2.0)*ellk(2.0*p))
      eta_2=(p^2.0/2.0)*(p^2.0+2.0*z[w]^2.0)
      eta_1=(2.0*!PI)^(-1.0)*(k1+2.0*eta_2*k0-0.25*(1.0+5.0*p^2.0+z[w]^2.0)*SQRT((1.0-a)*(b-1.0)))
      lambda_d=lambda4
      eta_d=eta_2
      IF (p LE z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
      ENDIF
      IF (p GT z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
      ENDIF
   ENDIF

   ;case VI
   IF (p EQ 0.5) and (z[w] EQ 0.5) THEN BEGIN
      ;print, (a-b)/a
      ;print, k^(-1.0)
      ;print, ellpic_bulirsch(((a-b)/a),k^(-1.0))
      lambda_d=(1.0/3.0)-(4.0/(9.0*!PI))
      eta_d=3.0/32.0
      IF (p LE z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
      ENDIF
      IF (p GT z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
      ENDIF
   ENDIF

   ;case VII
   IF (p GT 0.5) and (z[w] EQ p) THEN BEGIN
      ;print, (a-b)/a
      ;print, k^(-1.0)
      ;print, ellpic_bulirsch(((a-b)/a),k^(-1.0))
      lambda3=(1.0/3.0)+((16.0*p)/(9.0*!PI))*(2.0*p^2.0-1.0)*ellec(1.0/(2.0*p))-(((1.0-4.0*p^2.0)*(3.0-8.0*p^2.0))/(9.0*!PI*p))*ellk(1.0/(2.0*p))
      eta_2=(p^2.0/2.0)*(p^2.0+2.0*z[w]^2.0)
      eta_1=(2.0*!PI)^(-1.0)*(k1+2.0*eta_2*k0-0.25*(1.0+5.0*p^2.0+z[w]^2.0)*SQRT((1.0-a)*(b-1.0)))
      lambda_d=lambda3
      eta_d=eta_1
      IF (p LE z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
      ENDIF
      IF (p GT z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
      ENDIF
   ENDIF

   ;case VIII
   IF (p GT 0.5) and (z[w] GE ABS(1.0-p)) and (z[w] LT p) THEN BEGIN
      ;print, (a-b)/a
      ;print, k^(-1.0)
      ;print, ellpic_bulirsch(((a-b)/a),k^(-1.0))
      lambda1=(1.0/(9.0*!PI*SQRT(p*z[w])))*(((1.0-b)*(2.0*b+a-3.0)-3.0*q*(b-2.0))*ellk(k)+4.0*p*z[w]*(z[w]^2.0+7.0*p^2.0-4.0)*ellec(k)-(3.0*q/a)*ellpic_bulirsch(ABS((a-1.0)/a),k))
      eta_2=(p^2.0/2.0)*(p^2.0+2.0*z[w]^2.0)
      eta_1=(2.0*!PI)^(-1.0)*(k1+2.0*eta_2*k0-0.25*(1.0+5.0*p^2.0+z[w]^2.0)*SQRT((1.0-a)*(b-1.0)))
      lambda_d=lambda1
      eta_d=eta_1
      IF (p LE z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
      ENDIF
      IF (p GT z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
      ENDIF
   ENDIF

   ;case IX
   IF (p GT 0.0) and (p LT 1.0) and (z[w] GT 0.0) and (z[w] LT 0.5-ABS(p-0.5)) THEN BEGIN
      ;print, (a-b)/a
      ;print, k^(-1.0)
      ;print, ellpic_bulirsch(((a-b)/a),k^(-1.0))
      lambda2=(2.0/(9.0*!PI*SQRT(1.0-a)))*((1.0-5.0*z[w]^2.0+p^2.0+q^2.0)*ellk(k^(-1.0))+(1.0-a)*(z[w]^2.0+7.0*p^2.0-4.0)*ellec(k^(-1.0))-(3.0*q/a)* ellpic_bulirsch(ABS((a-b)/a),k^(-1.0)))
      eta_2=(p^2.0/2.0)*(p^2.0+2.0*z[w]^2.0)
      eta_1=(2.0*!PI)^(-1.0)*(k1+2.0*eta_2*k0-0.25*(1.0+5.0*p^2.0+z[w]^2.0)*SQRT((1.0-a)*(b-1.0)))
      lambda_d=lambda2
      eta_d=eta_2
      IF (p LE z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
      ENDIF
      IF (p GT z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
      ENDIF
   ENDIF

   ;case X
   IF (p GT 0.0) and (p LT 1.0) and (z[w] EQ 0.0) THEN BEGIN
      ;print, (a-b)/a
      ;print, k^(-1.0)
      ;print, ellpic_bulirsch(((a-b)/a),k^(-1.0))
      lambda6=-(2.0/3.0)*(1.0-p^2.0)^(3.0/2.0)
      eta_2=(p^2.0/2.0)*(p^2.0+2.0*z[w]^2.0)
      eta_1=(2.0*!PI)^(-1.0)*(k1+2.0*eta_2*k0-0.25*(1.0+5.0*p^2.0+z[w]^2.0)*SQRT((1.0-a)*(b-1.0)))
      lambda_d=lambda6
      eta_d=eta_2
      IF (p LE z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
      ENDIF
      IF (p GT z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
      ENDIF
   ENDIF

   ;case XI
   IF (p GT 1.0) and (z[w] GE 0.0) and (z[w] LT p-1.0) THEN BEGIN
      ;print, (a-b)/a
      ;print, k^(-1.0)
      ;print, ellpic_bulirsch(((a-b)/a),k^(-1.0))
      lambda_d=0.0
      eta_d=0.5
      IF (p LE z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d-c4*eta_d))
      ENDIF
      IF (p GT z[w]) THEN BEGIN
         F[w]=1.0-(4.0*Omega)^(-1.0)*((1.0-c2)*lambda_e+c2*(lambda_d+(2.0/3.0)-c4*eta_d))
      ENDIF
   ENDIF

ENDFOR


;convers�o de x para tempo
ttrans=periodo/(!PI*adivR)
meio = (Nx/2.0)-0.5
FOR u=0, Nx-1 DO BEGIN
   IF (u LT meio) THEN BEGIN
      curva_simulada[0,u]=(-1.0)*xs[u]*ttrans/2.0
   ENDIF
   IF (u GE meio) THEN BEGIN
      curva_simulada[0,u]=xs[u]*ttrans/2.0
   ENDIF
ENDFOR
curva_simulada[1,*]=F[*]

;reamostragem
curva_simulada_reamostrada[0,*]=curva_observada_eclipse[0,*]
curva_simulada_reamostrada[1,*]=INTERPOL(curva_simulada[1,*],curva_simulada[0,*],curva_simulada_reamostrada[0,*])

OPENW,1,curva_saida
FOR w=0, Nobs-1 DO BEGIN
   PRINTF,1,curva_simulada_reamostrada[*,w]
ENDFOR
close,1







END