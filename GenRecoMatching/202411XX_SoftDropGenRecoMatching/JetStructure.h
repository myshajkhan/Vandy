#ifndef JETSTRUCTURE_H   // match the file name, all caps, underscores
#define JETSTRUCTURE_H

// We are creating construction fuction "Jet".
struct Jet
{
    
    //Declaring the member varibales for Jet
    double PT;
    double Eta;
    double Phi;
    double E;
    double ZG;
    double RG;
    
    Jet(double pt, double eta, double phi, double e, double zg, double rg)
    : PT(pt), Eta(eta), Phi(phi), E(e), ZG(zg), RG(rg)
    {
        
    }
};

//Metric fuction is finding the angle between gen and reco jet
double Metric(Jet Gen, Jet Reco)
{
    double gen_px =  Gen.PT * cos(Gen.Phi); // Calculate the x-component of the gen jet's momentum using transverse momentum and azimuthal angle
    double gen_py =  Gen.PT * sin(Gen.Phi);  // Calculate the y-component of the gen jet's momentum
    double gen_pz =  Gen.PT * sinh(Gen.Eta); // Calculate the z-component of the gen jet's momentum using pseudorapidity
    double gen_p = Gen.PT * cosh(Gen.Eta); // Calculate the magnitude of the gen jet's momentum vector
    double reco_px =  Reco.PT * cos(Reco.Phi); // Calculate the x-component of the reco jet's momentum
    double reco_py =  Reco.PT * sin(Reco.Phi); // Calculate the y-component of the reco jet's momentum
    double reco_pz =  Reco.PT * sinh(Reco.Eta);  // Calculate the z-component of the reco jet's momentum
    double reco_p = Reco.PT * cosh(Reco.Eta); // Calculate the magnitude of the reco jet's momentum vector
    double gen_reco_dot = (gen_px * reco_px) + (gen_py * reco_py) +  (gen_pz * reco_pz) ;  // Compute the dot product between the gen and reco momentum vectors
    double cos_angle = gen_reco_dot / (gen_p *  reco_p) ;  // Compute the cosine of the angle between the two vectors using the dot product formula
    
    // avoing cos value larger than 1, so we dont get arc cos = nan
    if ( cos_angle > 0.9999999999 ) cos_angle = 0.9999999999 ;
    if ( cos_angle < -0.9999999999 ) cos_angle = -0.9999999999 ;
    return acos(cos_angle);
}
#endif // JETSTRUCTURE_H
