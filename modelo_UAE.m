function dCadt = modelo_UAE(t, Ca, e, ka_s, ka_f, a, K)
  % Ca(1) = Ca_s: concentración en el solvente
  % Ca(2) = Ca_f: concentración en la fruta
  Ca_s = Ca(1);
  Ca_f = Ca(2);

  % --- Condición de equilibrio en la interfaz ---
  %   Ca_s,i = K * Ca_f,i
  % --- Flujo continuo en la interfaz ---
  %   ka_s*(Ca_s,i - Ca_s) = -ka_f*(Ca_f,i - Ca_f)
  %
  % Sustituimos Ca_s,i = K*Ca_f,i en la segunda:
  %   ka_s*(K*Ca_f_i - Ca_s) = -ka_f*(Ca_f_i - Ca_f)
  % => resolvemos para Ca_f_i:
  Ca_f_i = (ka_s * Ca_s + ka_f * Ca_f) / (ka_s * K + ka_f);
  Ca_s_i = K * Ca_f_i;

  % --- Ecuaciones macroscópicas de transferencia (7) y (8) ---
  % (7) ? dCa_s/dt = ka_s * S * (Ca_s,i - Ca_s)
  % (8) (1-?) dCa_f/dt = ka_f * S * (Ca_f_i - Ca_f)
  dCa_s = (ka_s * a * (Ca_s_i - Ca_s)) / e;
  dCa_f = (ka_f *a * (Ca_f_i - Ca_f)) / (1 - e);

  dCadt = [dCa_s;dCa_f];
end
