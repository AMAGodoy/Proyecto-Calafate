function dCadt = modelo_antocianina(t, Ca, kc, K, a)
  Ca_f = Ca(1);  % Concentración en fruta
  Ca_s = Ca(2);  % Concentración en solvente

  dCa_f = -kc * a * (K*Ca_f - Ca_s);
  dCa_s =  kc * a * (K*Ca_f - Ca_s);

  dCadt = [dCa_f; dCa_s];
end
