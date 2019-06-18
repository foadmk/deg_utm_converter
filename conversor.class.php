<?php

class conversor {

    function __construct()
    {
        $this->a = 6378160.0 ;
        $this->f = 298.257223563;  //change to 298.25 for SAD69
        $this->k0 = 0.9996;
        $this->e2 = 2*(1/$this->f)-pow(1/$this->f,2);
        $this->e2linha = $this->e2/(1-$this->e2);
        $this->e1 = (1-sqrt(1-$this->e2))/(1+sqrt(1-$this->e2));
        $this->utm_z = 22;        
    }

    function utm_to_deg($utm_x, $utm_y)
    {

        $x = $utm_x - 500000;
        $y = $utm_y;

        $longOri = ($this->utm_z - 1 ) * 6 - 180 + 3;

        $M = $y/$this->k0;
        $mu = $M / ($this->a*(1-$this->e2/4-3*pow($this->e2,2)/64-5*pow($this->e2,3)/256));

        //echo "M: {$M}  /  mu: {$mu} <br>";

        $p1 = 3*$this->e1/2-27*pow($this->e1,3)/27;
        $p2 = 21*pow($this->e1,2)/16 + 55*pow($this->e1,4)/32;
        $p3 = 151*pow($this->e1,3)/96;

        $phiRad = $mu + sin(2*$mu)*$p1 + sin(4*$mu)*$p2 + sin(6*$mu)*$p3;
        $phi = 180.0 * $phiRad / pi();
        //echo "Phi Rad: {$phiRad}  /  Phi: {$phi} <br>";

        $N1 = $this->a / sqrt(1 - $this->e2*(pow(sin($phiRad),2)));
        $T1 = pow(tan($phiRad),2);
        $C1 = $this->e2linha*(pow(cos($phiRad),2));
        $R1 = $this->a*(1-$this->e2) / ( pow((1-$this->e2*pow(sin($phiRad),2)), 1.5) );
        $D = $x/($N1*$this->k0);

        //echo "N1: {$N1}  /  T1: {$T1} /  C1: {$C1} /  R1: {$R1} /  D: {$D} <br>";
        $latRad = $phiRad - ($N1*tan($phiRad)/$R1) * ( (pow($D,2)/2) - (5+3*$T1+10*$C1-4*pow($C1,2)-9*$this->e2linha)*(pow($D,4)/24) + (61+90*$T1+298*$C1+45*pow($T1,2)-252*$this->e2linha-3*pow($C1,2))*(pow($D,6)/720));
        $lonRad = ($D-(1+2*$T1+$C1)*pow($D,3)/6 + (5-2*$C1+28*$T1-3*pow($C1,2)+8*$this->e2linha+24*pow($T1,2))*pow($D,5)/120)/COS($phiRad);

        $lat = $latRad * 180 / pi();
        $long = $lonRad * 180 / pi() + $longOri;

        return ["lat" => $lat, "long" => $long];

    }

    function deg_to_utm($lat, $long) {

        $latRad = $lat * pi() / 180.0;
        $longRad = $long * pi() / 180.0;

        $longOri = ($this->utm_z - 1 ) * 6 - 180 + 3;
        $longOriRad = $longOri * pi() / 180.0;

        $N1 = $this->a / sqrt(1 - $this->e2*(pow(sin($latRad),2)));
        $T1 = pow(tan($latRad),2);
        $C1 = $this->e2linha*(pow(cos($latRad),2));

        $P1 = ($longRad-$longOriRad)*cos($latRad);
        $P2 = $this->a*((1-$this->e2/4-3*(pow($this->e2,2))/64-5*(pow($this->e2,3))/256)*$latRad - (3*$this->e2/8+3*(pow($this->e2,2))/32+45*(pow($this->e2,3))/1024)*sin(2*$latRad) + (15*(pow($this->e2,2))/256+45*(pow($this->e2,3))/1024)*sin(4*$latRad) - (35*(pow($this->e2,3))/3072)*sin(6*$latRad));

        //echo "N1: {$N1}  /  T1: {$T1} /  C1: {$C1} /  P1: {$P1} /  P2: {$P2} <br>";

        $utm_x = $this->k0*$N1*($P1 + ((1-$T1+$C1)*pow($P1,3))/6 + ((5-18*$T1+pow($T1,2)+72*$C1-58*$this->e2linha)*pow($P1,5))/120)+500000;

        $utm_y = $this->k0*($P2+$N1*tan($latRad)*(pow($P1,2)/2 + ((5-$T1+9*$C1+4*pow($C1,2))*pow($P1,4))/24 + ((61-58*$T1+pow($T1,2)+600*$C1-330*$this->e2linha)*pow($P1,6))/720));

        //echo "UTM_X: {$utm_x}  /  UTM_Y: {$utm_y} <br>";

        return ["utm_x" => $utm_x, "utm_y" => $utm_y];

    }
}

//Testing
$conv = new conversor();
print_r($deg = $conv->utm_to_deg(155244.82, 7154038.54));
print_r($utm = $conv->deg_to_utm($deg['lat'], $deg['long']));

?>
