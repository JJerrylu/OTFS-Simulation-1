% =========================================================================
% Book title: OTFS Modulation: Theory and Applications
% Author: Saif Khan Mohammed, Ronny Hadani and Ananthanarayanan
% Chockalingam
% Editor: Wiley-IEEE Press
% Date: Oct. 2024
% Description: This Matlab program computes value of the Root raised
% cosine pulse rrc_{\nu}(\nu) at \nu and returns it.
% See equation (2.82) in the book for the expression of the RRC pulse
% =========================================================================

function nuval = wnu(nu, beta_nu, constbetanu, abstolnu)

if (( abs(abs(beta_nu*4*nu) - 1) < abstolnu))
    
    nuval = constbetanu;
else
    nuval =  ( (4*beta_nu/pi)*cos((1 + beta_nu)*pi*nu) +  (1 - beta_nu)*sinc((1 - beta_nu)*nu)) ./ (1 - (4*beta_nu*nu).*(4*beta_nu*nu))   ;
end



end


