function [ x, y, z, w ] = ld1730 ( )

%*****************************************************************************80
%
%% ld1730() computes the 1730 point Lebedev angular grid.
%
%  Modified:
%
%    14 September 2010
%
%  Author:
%
%    Dmitri Laikov
%
%  Reference:
%
%    Vyacheslav Lebedev, Dmitri Laikov,
%    A quadrature formula for the sphere of the 131st
%    algebraic order of accuracy,
%    Russian Academy of Sciences Doklady Mathematics,
%    Volume 59, Number 3, 1999, pages 477-481.
%
%  Output:
%
%    real X(N), Y(N), Z(N), W(N), the coordinates
%    and weights of the points.
%
  n = 0;
  x = zeros(1730,1);
  y = zeros(1730,1);
  z = zeros(1730,1);
  w = zeros(1730,1);
  a = 0.0;
  b = 0.0;
  v = 0.6309049437420976E-04;
  [ n, x, y, z, w ] = gen_oh ( 1, n, a, b, v, x, y, z, w );
  v = 0.6398287705571748E-03;
  [ n, x, y, z, w ] = gen_oh ( 2, n, a, b, v, x, y, z, w );
  v = 0.6357185073530720E-03;
  [ n, x, y, z, w ] = gen_oh ( 3, n, a, b, v, x, y, z, w );
  a = 0.2860923126194662E-01;
  v = 0.2221207162188168E-03;
  [ n, x, y, z, w ] = gen_oh ( 4, n, a, b, v, x, y, z, w );
  a = 0.7142556767711522E-01;
  v = 0.3475784022286848E-03;
  [ n, x, y, z, w ] = gen_oh ( 4, n, a, b, v, x, y, z, w );
  a = 0.1209199540995559;
  v = 0.4350742443589804E-03;
  [ n, x, y, z, w ] = gen_oh ( 4, n, a, b, v, x, y, z, w );
  a = 0.1738673106594379;
  v = 0.4978569136522127E-03;
  [ n, x, y, z, w ] = gen_oh ( 4, n, a, b, v, x, y, z, w );
  a = 0.2284645438467734;
  v = 0.5435036221998053E-03;
  [ n, x, y, z, w ] = gen_oh ( 4, n, a, b, v, x, y, z, w );
  a = 0.2834807671701512;
  v = 0.5765913388219542E-03;
  [ n, x, y, z, w ] = gen_oh ( 4, n, a, b, v, x, y, z, w );
  a = 0.3379680145467339;
  v = 0.6001200359226003E-03;
  [ n, x, y, z, w ] = gen_oh ( 4, n, a, b, v, x, y, z, w );
  a = 0.3911355454819537;
  v = 0.6162178172717512E-03;
  [ n, x, y, z, w ] = gen_oh ( 4, n, a, b, v, x, y, z, w );
  a = 0.4422860353001403;
  v = 0.6265218152438485E-03;
  [ n, x, y, z, w ] = gen_oh ( 4, n, a, b, v, x, y, z, w );
  a = 0.4907781568726057;
  v = 0.6323987160974212E-03;
  [ n, x, y, z, w ] = gen_oh ( 4, n, a, b, v, x, y, z, w );
  a = 0.5360006153211468;
  v = 0.6350767851540569E-03;
  [ n, x, y, z, w ] = gen_oh ( 4, n, a, b, v, x, y, z, w );
  a = 0.6142105973596603;
  v = 0.6354362775297107E-03;
  [ n, x, y, z, w ] = gen_oh ( 4, n, a, b, v, x, y, z, w );
  a = 0.6459300387977504;
  v = 0.6352302462706235E-03;
  [ n, x, y, z, w ] = gen_oh ( 4, n, a, b, v, x, y, z, w );
  a = 0.6718056125089225;
  v = 0.6358117881417972E-03;
  [ n, x, y, z, w ] = gen_oh ( 4, n, a, b, v, x, y, z, w );
  a = 0.6910888533186254;
  v = 0.6373101590310117E-03;
  [ n, x, y, z, w ] = gen_oh ( 4, n, a, b, v, x, y, z, w );
  a = 0.7030467416823252;
  v = 0.6390428961368665E-03;
  [ n, x, y, z, w ] = gen_oh ( 4, n, a, b, v, x, y, z, w );
  a = 0.8354951166354646E-01;
  v = 0.3186913449946576E-03;
  [ n, x, y, z, w ] = gen_oh ( 5, n, a, b, v, x, y, z, w );
  a = 0.2050143009099486;
  v = 0.4678028558591711E-03;
  [ n, x, y, z, w ] = gen_oh ( 5, n, a, b, v, x, y, z, w );
  a = 0.3370208290706637;
  v = 0.5538829697598626E-03;
  [ n, x, y, z, w ] = gen_oh ( 5, n, a, b, v, x, y, z, w );
  a = 0.4689051484233963;
  v = 0.6044475907190476E-03;
  [ n, x, y, z, w ] = gen_oh ( 5, n, a, b, v, x, y, z, w );
  a = 0.5939400424557334;
  v = 0.6313575103509012E-03;
  [ n, x, y, z, w ] = gen_oh ( 5, n, a, b, v, x, y, z, w );
  a = 0.1394983311832261;
  b = 0.4097581162050343E-01;
  v = 0.4078626431855630E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.1967999180485014;
  b = 0.8851987391293348E-01;
  v = 0.4759933057812725E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.2546183732548967;
  b = 0.1397680182969819;
  v = 0.5268151186413440E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.3121281074713875;
  b = 0.1929452542226526;
  v = 0.5643048560507316E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.3685981078502492;
  b = 0.2467898337061562;
  v = 0.5914501076613073E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.4233760321547856;
  b = 0.3003104124785409;
  v = 0.6104561257874195E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.4758671236059246;
  b = 0.3526684328175033;
  v = 0.6230252860707806E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.5255178579796463;
  b = 0.4031134861145713;
  v = 0.6305618761760796E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.5718025633734589;
  b = 0.4509426448342351;
  v = 0.6343092767597889E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.2686927772723415;
  b = 0.4711322502423248E-01;
  v = 0.5176268945737826E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.3306006819904809;
  b = 0.9784487303942695E-01;
  v = 0.5564840313313692E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.3904906850594983;
  b = 0.1505395810025273;
  v = 0.5856426671038980E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.4479957951904390;
  b = 0.2039728156296050;
  v = 0.6066386925777091E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.5027076848919780;
  b = 0.2571529941121107;
  v = 0.6208824962234458E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.5542087392260217;
  b = 0.3092191375815670;
  v = 0.6296314297822907E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.6020850887375187;
  b = 0.3593807506130276;
  v = 0.6340423756791859E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.4019851409179594;
  b = 0.5063389934378671E-01;
  v = 0.5829627677107342E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.4635614567449800;
  b = 0.1032422269160612;
  v = 0.6048693376081110E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.5215860931591575;
  b = 0.1566322094006254;
  v = 0.6202362317732461E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.5758202499099271;
  b = 0.2098082827491099;
  v = 0.6299005328403779E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.6259893683876795;
  b = 0.2618824114553391;
  v = 0.6347722390609353E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.5313795124811891;
  b = 0.5263245019338556E-01;
  v = 0.6203778981238834E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.5893317955931995;
  b = 0.1061059730982005;
  v = 0.6308414671239979E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.6426246321215801;
  b = 0.1594171564034221;
  v = 0.6362706466959498E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  a = 0.6511904367376113;
  b = 0.5354789536565540E-01;
  v = 0.6375414170333233E-03;
  [ n, x, y, z, w ] = gen_oh ( 6, n, a, b, v, x, y, z, w );
  
  return
end