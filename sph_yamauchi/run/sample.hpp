#pragma once

#include "hdr_heos.hpp"
#include "hdr_util.hpp"

static PS::U64 convertF64ToU64(PS::F64 val) {
    union converter {
        PS::F64 f;
        PS::U64 u;
    };
    union converter var;
        var.f = val;
        return var.u;
}
static PS::F64 convertU64ToF64(PS::U64 val) {
    union converter {
        PS::F64 f;
        PS::U64 u;
    };
    union converter var;
    var.u = val;
    return var.f;
}

class HelmholtzGas : public SPH {
public:
    PS::S64 istar;
    PS::F64 temp;
    PS::F64 entr;
    PS::F64 dnuc;
    PS::F64 enuc;
    NR::Nucleon cmps;
    PS::F64 pot3;
    PS::F64 gcor3, dgcor3;
    PS::F64 tempmax[3];
    bool flagwrite;

    static const PS::F64 nidx;
    static const PS::F64 ninv;
    static const PS::F64 kcon;
    static const PS::F64 kconInThisUnit;

    HelmholtzGas() {
        id     = 0;
        mass   = 0.;
        pos    = 0.;
        vel    = 0.;
        vel2   = 0.;
        acc    = 0.;
        uene   = 0.;
        uene2  = 0.;
        udot   = 0.;
        alph   = 0.;
        alph2  = 0.;
        adot   = 0.;
        alphu  = 0.;
        alphu2 = 0.;    
        adotu  = 0.;
        diffu  = 0.;
        vol    = 0.;
        ksr    = 0.;
        rs     = 0.;
        dens   = 0.;
        pres   = 0.;
        vsnd   = 0.;
        divv   = 0.;
        rotv   = 0.;
        bswt   = 0.;
        np     = 0;
        grdh   = 0.;
        vsmx   = 0.;
        pot    = 0.;
        acch   = 0.;
        accg1  = 0.;
        accg2  = 0.;
        eta    = 0.;
        abar   = 0.;
        zbar   = 0.;
        umin   = 0.;
        istar  = 0;
        temp   = 0.;
        entr   = 0.;
        dnuc   = 0.;
        enuc   = 0.;
        pot3   = 0.;
        gcor3  = 0.;
        dgcor3 = 0.;
        flagwrite = false;
    }

    void readAscii(FILE * fp) {
        using namespace CodeUnit;
        fscanf(fp, "%lld%lld%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
               &this->id, &this->istar, &this->mass,        // 3
               &this->pos[0], &this->pos[1], &this->pos[2], // 6
               &this->vel[0], &this->vel[1], &this->vel[2], // 9
               &this->uene,   &this->alph,   &this->alphu,  // 12
               &this->ksr);                                 // 13
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {  // 26
            fscanf(fp, "%lf", &this->cmps[k]);
        }
        this->mass *= UnitOfMassInv;
        this->pos  *= UnitOfLengthInv;
        this->vel  *= UnitOfVelocityInv;
        this->uene *= UnitOfEnergyInv;
        this->ksr  *= UnitOfLengthInv;
    }
    
    void writeAscii(FILE *fp) const {
        using namespace CodeUnit;
        PS::F64    tmass = this->mass * UnitOfMass;
        PS::F64vec tpos  = this->pos  * UnitOfLength;
        PS::F64vec tvel  = this->vel  * UnitOfVelocity;
        PS::F64vec tacc  = this->acc  * UnitOfAcceleration;
        PS::F64    tuene = this->uene * UnitOfEnergy;
        PS::F64    tksr  = this->ksr  * UnitOfLength;
        PS::F64    tdens = this->dens * UnitOfDensity;
        PS::F64    tvsnd = this->vsnd * UnitOfVelocity;
        PS::F64    tpres = this->pres * UnitOfPressure;
        PS::F64    tdivv = this->divv * UnitOfTimeInv; // divv [s^-1]
        PS::F64    trotv = this->rotv * UnitOfTimeInv; // rotv [s^-1]
        PS::F64    tpot  = this->pot  * UnitOfEnergy;
        PS::F64    tenuc = this->enuc * UnitOfEnergy;
        PS::F64    tvsmx = this->vsmx * UnitOfVelocity;
        PS::F64    tudot = this->udot * UnitOfEnergy * UnitOfTimeInv;
        PS::F64    tdnuc = this->dnuc * UnitOfEnergy;
        PS::F64    tentr = this->entr * UnitOfEnergy;
        PS::F64    tpot3 = this->pot3 * UnitOfEnergy;
        fprintf(fp, "%6d %2d %+e", this->id, this->istar, tmass);    //  3
        fprintf(fp, " %+e %+e %+e", tpos[0], tpos[1], tpos[2]);      //  6
        fprintf(fp, " %+e %+e %+e", tvel[0], tvel[1], tvel[2]);      //  9
        fprintf(fp, " %+e %+e %+e", tacc[0], tacc[1], tacc[2]);      // 12
        fprintf(fp, " %+e %+e %+e", tuene, this->alph, this->alphu); // 15
        fprintf(fp, " %+e %+e %6d", tdens, tksr,  this->np);         // 18
        fprintf(fp, " %+e %+e %+e", tvsnd, tpres, this->temp);       // 21
        fprintf(fp, " %+e %+e %+e", tdivv, trotv, this->bswt);       // 24
        fprintf(fp, " %+e %+e %+e", tpot, this->abar, this->zbar);   // 27
        fprintf(fp, " %+e",         tenuc);                          // 28
        fprintf(fp, " %+e %+e %+e", tvsmx, tudot, tdnuc);            // 31
        for(PS::S32 k = 0; k < NuclearReaction::NumberOfNucleon; k++) { // 32 -- 44
            fprintf(fp, " %+.3e", this->cmps[k]);
        }
        if(RP::FlagDamping == 1) {
            PS::F64vec tomg = RP::RotationalVelocity * UnitOfTimeInv;
            PS::F64vec tvec = tomg ^ tpos;
            fprintf(fp, " %+e", tpot - 0.5 * (tvec * tvec));         // 45
        }
        fprintf(fp, " %+e", tpot3);                                  // 45
        fprintf(fp, " %+e %+e %+e", this->tempmax[0], this->tempmax[1] * UnitOfDensity,
                this->tempmax[2]);                                   // 46 -- 48
        fprintf(fp, " %+e", tentr);                                  // 49
        fprintf(fp, "\n");
    }

    void readHexa(FILE *fp) {
        PS::F64 (*cvt)(PS::U64) = convertU64ToF64;
        PS::S32 s0, s1;
        PS::U64 umass;
        PS::U64 upos[3], uvel[3], uacc[3], uacch[3], uaccg1[3], uaccg2[3];
        PS::U64 uuene, uudot;
        PS::U64 ualph, uadot, ualphu, uadotu;
        PS::U64 udens, upres, uvsnd, utemp;
        PS::U64 udivv, urotv, ubswt;
        PS::U64 uksr, ugrdh, uvsmx;
        PS::U64 upot, uenuc, udnuc;
        PS::U64 uabar, uzbar, ugcor3;
        PS::U64 ucmps[NR::NumberOfNucleon];
        PS::U64 utempmax[3];

        fscanf(fp, "%lld %lld %llx", &this->id, &this->istar, &umass);
        fscanf(fp, "%llx %llx %llx", &upos[0], &upos[1], &upos[2]);
        fscanf(fp, "%llx %llx %llx", &uvel[0], &uvel[1], &uvel[2]);
        fscanf(fp, "%llx %llx %llx", &uacc[0], &uacc[1], &uacc[2]);
        fscanf(fp, "%llx %llx %llx", &uacch[0], &uacch[1], &uacch[2]);
        fscanf(fp, "%llx %llx %llx", &uaccg1[0], &uaccg1[1], &uaccg1[2]);
        fscanf(fp, "%llx %llx %llx", &uaccg2[0], &uaccg2[1], &uaccg2[2]);
        fscanf(fp, "%llx %llx", &uuene, &uudot);
        fscanf(fp, "%llx %llx", &ualph, &uadot);
        fscanf(fp, "%llx %llx", &ualphu, &uadotu);
        fscanf(fp, "%llx %llx %llx %llx", &udens, &upres, &uvsnd, &utemp);
        fscanf(fp, "%llx %llx %llx", &udivv, &urotv, &ubswt);
        fscanf(fp, "%llx %llx %llx", &uksr,  &ugrdh, &uvsmx);
        fscanf(fp, "%llx %llx %llx", &upot, &uenuc, &udnuc);
        fscanf(fp, "%llx %llx %llx", &uabar, &uzbar, &ugcor3);
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            fscanf(fp, "%llx", &ucmps[k]);
        }
        fscanf(fp, "%llx %llx %llx", &utempmax[0], &utempmax[1], &utempmax[2]);

        this->mass = cvt(umass);
        this->pos[0] = cvt(upos[0]);
        this->pos[1] = cvt(upos[1]);
        this->pos[2] = cvt(upos[2]);
        this->vel[0] = cvt(uvel[0]);
        this->vel[1] = cvt(uvel[1]);
        this->vel[2] = cvt(uvel[2]);
        this->acc[0] = cvt(uacc[0]);
        this->acc[1] = cvt(uacc[1]);
        this->acc[2] = cvt(uacc[2]);
        this->acch[0] = cvt(uacch[0]);
        this->acch[1] = cvt(uacch[1]);
        this->acch[2] = cvt(uacch[2]);
        this->accg1[0] = cvt(uaccg1[0]);
        this->accg1[1] = cvt(uaccg1[1]);
        this->accg1[2] = cvt(uaccg1[2]);
        this->accg2[0] = cvt(uaccg2[0]);
        this->accg2[1] = cvt(uaccg2[1]);
        this->accg2[2] = cvt(uaccg2[2]);
        this->uene  = cvt(uuene);
        this->udot  = cvt(uudot);
        this->alph  = cvt(ualph);
        this->adot  = cvt(uadot);
        this->alphu = cvt(ualphu);
        this->adotu = cvt(uadotu);
        this->dens  = cvt(udens);
        this->pres  = cvt(upres);
        this->vsnd  = cvt(uvsnd);
        this->temp  = cvt(utemp);
        this->divv  = cvt(udivv);
        this->rotv  = cvt(urotv);
        this->bswt  = cvt(ubswt);
        this->ksr   = cvt(uksr);
        this->grdh  = cvt(ugrdh);
        this->vsmx  = cvt(uvsmx);
        this->pot   = cvt(upot);        
        this->enuc  = cvt(uenuc);        
        this->dnuc  = cvt(udnuc);        
        this->abar  = cvt(uabar);
        this->zbar  = cvt(uzbar);
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            this->cmps[k] = cvt(ucmps[k]);
        }
        this->gcor3      = cvt(ugcor3);
        this->tempmax[0] = cvt(utempmax[0]);
        this->tempmax[1] = cvt(utempmax[1]);
        this->tempmax[2] = cvt(utempmax[2]);
    }

    void writeHexa(FILE *fp) const {
        PS::U64 (*cvt)(PS::F64) = convertF64ToU64;
        fprintf(fp, "%6lld %2lld %llx", this->id, this->istar, cvt(this->mass));
        fprintf(fp, " %llx %llx %llx", cvt(this->pos[0]), cvt(this->pos[1]), cvt(this->pos[2]));
        fprintf(fp, " %llx %llx %llx", cvt(this->vel[0]), cvt(this->vel[1]), cvt(this->vel[2]));
        fprintf(fp, " %llx %llx %llx", cvt(this->acc[0]), cvt(this->acc[1]), cvt(this->acc[2]));
        fprintf(fp, " %llx %llx %llx", cvt(this->acch[0]), cvt(this->acch[1]), cvt(this->acch[2]));
        fprintf(fp, " %llx %llx %llx",
                cvt(this->accg1[0]), cvt(this->accg1[1]), cvt(this->accg1[2]));
        fprintf(fp, " %llx %llx %llx",
                cvt(this->accg2[0]), cvt(this->accg2[1]), cvt(this->accg2[2]));
        fprintf(fp, " %llx %llx", cvt(this->uene), cvt(this->udot));
        fprintf(fp, " %llx %llx", cvt(this->alph), cvt(this->adot));
        fprintf(fp, " %llx %llx", cvt(this->alphu), cvt(this->adotu));
        fprintf(fp, " %llx %llx", cvt(this->dens), cvt(this->pres));
        fprintf(fp, " %llx %llx", cvt(this->vsnd), cvt(this->temp));
        fprintf(fp, " %llx %llx %llx", cvt(this->divv), cvt(this->rotv), cvt(this->bswt));
        fprintf(fp, " %llx %llx %llx", cvt(this->ksr), cvt(this->grdh), cvt(this->vsmx));
        fprintf(fp, " %llx %llx %llx", cvt(this->pot), cvt(this->enuc), cvt(this->dnuc));
        fprintf(fp, " %llx %llx %llx", cvt(this->abar), cvt(this->zbar), cvt(this->gcor3));
        for(PS::S32 k = 0; k < NR::NumberOfNucleon; k++) {
            fprintf(fp, " %llx", cvt(this->cmps[k]));
        }
        fprintf(fp, " %llx %llx %llx", cvt(this->tempmax[0]), cvt(this->tempmax[1]),
                cvt(this->tempmax[2]));
        fprintf(fp, "\n");
    }
    
    void referEquationOfState() {
	this->pres = this->ninv * this->dens * this->uene;
	this->vsnd = sqrt((1. + this->ninv) * this->pres / this->dens);
	this->umin = 0.;
	if(RP::FlagDamping == 1) {
	    this->pres = this->kconInThisUnit * pow(this->dens, (1. + this->ninv));
	}
    }
    
    PS::F64 calcTimestep() {
        using namespace CodeUnit;
        PS::F64 tceff = RP::CoefficientOfTimestep;
        PS::F64 dth = tceff * this->ksr / (this->vsmx * SK::ksrh);
        PS::F64 dtu = RP::MaximumTimestep;
        PS::F64 dtc = (this->divv < 0.) ? (- tceff / this->divv) : RP::MaximumTimestep;
        PS::F64vec grv = this->accg1 + this->accg2;
        PS::F64 dtg = tceff * sqrt(this->ksr / sqrt(grv * grv));
        return std::min(std::min(dth, dtu), std::min(dtg, dtc));
    }

    inline void calcNewtonGravity(PS::F64vec & accg3,
                                  PS::F64 & pot3) {
        PS::F64    rinv  = 1. / sqrt(this->pos * this->pos);
        accg3 = (- CodeUnit::GravityConstantInThisUnit
                 * CodeUnit::BlackHoleMassInThisUnit * rinv * rinv * rinv) * this->pos;
        pot3  = - CodeUnit::GravityConstantInThisUnit * CodeUnit::BlackHoleMassInThisUnit * rinv;
    }

    inline void addAdditionalForce() {
        if(RP::FlagBinary == 2) {
            PS::F64vec accg3;
            PS::F64    pot3;
            PS::F64    dgcor3 = 0.0;
            if(RP::FlagPotential == 0) {
                calcNewtonGravity(accg3, pot3);
            } else {
                fprintf(stderr, "Unknown potential!\n");
                PS::Abort();
            }
            this->accg1 += accg3;
            this->acc   += accg3;
            this->pot   += (2. * pot3);
            this->pot3   = pot3;
            this->dgcor3 = dgcor3;
        }
    }

    inline void addAdditionalForceDamping1() {
        this->acc  -= RP::RotationalVelocity ^ (RP::RotationalVelocity ^ this->pos)
            + 2. * (RP::RotationalVelocity ^ this->vel);
    }

    PS::F64 calcEnergy() {
        return this->mass * (0.5 * this->vel * this->vel + this->uene + 0.5 * this->pot
                             - this->gcor3);
    }

    PS::F64 calcEnergyDamping1() {
	PS::F64vec tvec = RP::RotationalVelocity ^ this->pos;
        return this->mass * (0.5 * this->vel * this->vel + this->uene
                             + 0.5 * (this->pot - tvec * tvec));
    }

    PS::F64vec calcMomentum() {
        return this->mass * this->vel;
    }

    void predict(PS::F64 dt) {
        this->pos    = this->pos   +       this->vel   * dt  + 0.5 * this->acc * dt * dt;
        this->vel2   = this->vel   + 0.5 * this->acc   * dt;
        this->vel    = this->vel   +       this->acc   * dt;
        this->uene2  = this->uene  + 0.5 * this->udot  * dt;
        this->uene   = this->uene  +       this->udot  * dt  +       this->dnuc;
        this->alph2  = this->alph  + 0.5 * this->adot  * dt;
        this->alph   = this->alph  +       this->adot  * dt;
        this->alphu2 = this->alphu + 0.5 * this->adotu * dt;
        this->alphu  = this->alphu +       this->adotu * dt;
        this->gcor3  = this->gcor3 + 0.5 * this->dgcor3 * dt;
    }

    void correct(PS::F64 dt) {
        this->vel   = this->vel2   + 0.5 * this->acc   * dt;
        this->uene  = this->uene2  + 0.5 * this->udot  * dt  +       this->dnuc;
        this->alph  = this->alph2  + 0.5 * this->adot  * dt;
        this->alphu = this->alphu2 + 0.5 * this->adotu * dt;
        this->gcor3 = this->gcor3 + 0.5 * this->dgcor3 * dt;
    }

    void correctDamping1(PS::F64 dt) {
        //this->acc  -= this->vel * 0.1;
	//this->acc  -= this->vel * 0.01;
	this->acc  -= this->vel * 0.001;
        this->vel   = this->vel2   + 0.5 * this->acc  * dt;
        this->uene  = this->uene2  + 0.5 * this->udot * dt;
        this->alph  = this->alph2  + 0.5 * this->adot * dt;
        this->alphu = this->alphu2 + 0.5 * this->adotu * dt;
	this->uene  = this->pres / (this->ninv * this->dens);
    }

    void calcAlphaDot() {
        PS::F64 src    = std::max((- divv * (RP::AlphaMaximum - this->alph)), 0.);
        PS::F64 tauinv = (0.25 * SK::ksrh * this->vsnd) / this->ksr;
        this->adot     = - (this->alph - RP::AlphaMinimum) * tauinv + src;
        this->adotu    = 0.0;
	PS::F64 uene   = this->uene;
	PS::F64 srcu   = this->ksr * SK::ksrhinv * fabs(this->diffu)
            / sqrt(uene + RP::EpsilonOfInternalEnergy)
            * std::max(RP::AlphuMaximum - this->alphu, 0.);
        this->adotu    = - (this->alphu - RP::AlphuMinimum) * tauinv + srcu;
    }

};

const PS::F64 HelmholtzGas::nidx = 2.8;
const PS::F64 HelmholtzGas::ninv = 1. / HelmholtzGas::nidx;
const PS::F64 HelmholtzGas::kcon = 9.7028525e14;
const PS::F64 HelmholtzGas::kconInThisUnit = HelmholtzGas::kcon * CodeUnit::UnitOfPressureInv
    / pow(CodeUnit::UnitOfDensityInv, (1. + HelmholtzGas::ninv));

namespace RunParameter {    

    PS::S64 readAscii(const char * const filename) {
        using namespace CodeUnit;

        FILE *fp = fopen(filename, "r");
        if(fp == NULL) {
            fprintf(stderr, "Not found header file %s!\n", filename);
            PS::Abort();
        }
        fscanf(fp, "%lf%lf%lf%lf", &Time, &TimeEnd, &TimestepAscii, &TimestepHexa);
        fscanf(fp, "%lld%lld", &FlagDamping, &FlagNuclear);
        fscanf(fp, "%lf%lf", &AlphaMaximum, &AlphaMinimum);
        fscanf(fp, "%lf%lf", &AlphuMaximum, &AlphuMinimum);
        fscanf(fp, "%lf", &CoefficientOfTimestep);
        fscanf(fp, "%lld", &NumberOfParticle);
	if(FlagDamping == 1) {
	    fscanf(fp, "%lf", &RotationalVelocity[2]);
	}
        fclose(fp);

        Time          *= UnitOfTimeInv;
        TimeEnd       *= UnitOfTimeInv;
        TimestepAscii *= UnitOfTimeInv;
        TimestepHexa  *= UnitOfTimeInv;
    }


    template <class Tdinfo,
              class Tsph,
              class Tbhns,
              class Tmsls>
    void readHexa(FILE *fp,                         
                  Tdinfo & dinfo,
                  Tsph & sph,
                  Tbhns & bhns,
                  Tmsls & msls) {
        PS::F64 (*cvt)(PS::U64) = convertU64ToF64;        
        PS::U64 utime, utend, udtime, udta, udth;
        PS::U64 ualphamax, ualphamin, ualphumax, ualphumin;
        PS::U64 uksrmax, uepsu, uenuc, ueabs;
        PS::U64 urtime, urtimeinv, utcoeff;
        PS::U64 urv[3];
        PS::S64 nstep, nptcl, nmsls, nbhns, fdivide, fdamp, fnuc, fbin, fpot;
        PS::S64 nascii, nhexa;
        PS::U64 umbh, umbhunit;
        fscanf(fp, "%llx %llx %llx %llx %llx %lld", &utime, &utend, &udtime, &udta, &udth, &nstep);
        fscanf(fp, "%lld %lld", &nascii, &nhexa);
        fscanf(fp, "%llx %llx %llx %llx", &ualphamax, &ualphamin, &ualphumax, &ualphumin);
        fscanf(fp, "%llx %llx", &uksrmax, &uepsu);
        fscanf(fp, "%llx %llx", &urtime, &urtimeinv);
        fscanf(fp, "%lld %lld %lld %llx", &nptcl, &nmsls, &nbhns, &utcoeff);
        fscanf(fp, "%llx %llx %llx", &urv[0], &urv[1], &urv[2]);
        fscanf(fp, "%llx %llx %lld %lld %lld %lld %lld", &uenuc, &ueabs, &fdivide, &fdamp, &fnuc, &fbin, &fpot);
        fscanf(fp, "%llx %llx", &umbh, &umbhunit);
        Time          = cvt(utime);
        TimeEnd       = cvt(utend);
        Timestep      = cvt(udtime);
        TimestepAscii = cvt(udta);
        TimestepHexa  = cvt(udth);
        NumberOfStep  = nstep;
        if(RP::Time - (PS::S64)(RP::Time / RP::TimestepAscii) * RP::TimestepAscii == 0.) {
            NumberOfAscii = nascii - 1;
        } else {
            NumberOfAscii = nascii;
        }
        NumberOfHexa  = nhexa  - 1;
        AlphaMaximum  = cvt(ualphamax);
        AlphaMinimum  = cvt(ualphamin);
        AlphuMaximum  = cvt(ualphumax);
        AlphuMinimum  = cvt(ualphumin);
        KernelSupportRadiusMaximum = cvt(uksrmax);
        EpsilonOfInternalEnergy    = cvt(uepsu);
        ReductionTime    = cvt(urtime);
        ReductionTimeInv = cvt(urtimeinv);
        sph.setNumberOfParticleLocal(nptcl);
        msls.setNumberOfParticleLocal(nmsls);
        bhns.setNumberOfParticleLocal(nbhns);
        CoefficientOfTimestep = cvt(utcoeff);
        RotationalVelocity[0] = cvt(urv[0]);
        RotationalVelocity[1] = cvt(urv[1]);
        RotationalVelocity[2] = cvt(urv[2]);
        NuclearEnergyTotal    = cvt(uenuc);
        AbsorbedEnergyTotal   = cvt(ueabs);
        FlagDivideFile        = fdivide;
        FlagDamping           = fdamp;
        FlagNuclear           = fnuc;
        FlagBinary            = fbin;
        FlagPotential         = fpot;
        CodeUnit::BlackHoleMass           = cvt(umbh);
        CodeUnit::BlackHoleMassInThisUnit = cvt(umbhunit);

        for(PS::S32 i = 0; i < PS::Comm::getNumberOfProc(); i++) {
            PS::F64ort domain;
            PS::U64    u0, u1, u2, u3, u4, u5;
            fscanf(fp, "%llx %llx %llx %llx %llx %llx\n", &u0, &u1, &u2, &u3, &u4, &u5);
            domain.low_[0]  = cvt(u0);
            domain.low_[1]  = cvt(u1);
            domain.low_[2]  = cvt(u2);
            domain.high_[0] = cvt(u3);
            domain.high_[1] = cvt(u4);
            domain.high_[2] = cvt(u5);
            dinfo.setPosDomain(i, domain);
        }
    }

    template <class Tdinfo,
              class Tsph,
              class Tbhns,
              class Tmsls>
    void writeHexa(FILE *fp,                          
                   Tdinfo & dinfo,
                   Tsph & sph,
                   Tbhns & bhns,
                   Tmsls & msls) {
        PS::U64 (*cvt)(PS::F64) = convertF64ToU64;
        fprintf(fp, "%llx %llx %llx %llx %llx %lld\n",
                cvt(Time), cvt(TimeEnd), cvt(Timestep),
                cvt(TimestepAscii), cvt(TimestepHexa), NumberOfStep);
        fprintf(fp, "%lld %lld\n", NumberOfAscii, NumberOfHexa);
        fprintf(fp, "%llx %llx %llx %llx\n",
                cvt(AlphaMaximum), cvt(AlphaMinimum), cvt(AlphuMaximum), cvt(AlphuMinimum));
        fprintf(fp, "%llx %llx %llx %llx\n",
                cvt(KernelSupportRadiusMaximum), cvt(EpsilonOfInternalEnergy),
                cvt(ReductionTime), cvt(ReductionTimeInv));
        fprintf(fp, "%lld %lld %lld %llx\n", sph.getNumberOfParticleLocal(),
                msls.getNumberOfParticleLocal(), bhns.getNumberOfParticleLocal(),
                cvt(CoefficientOfTimestep));
        fprintf(fp, "%llx %llx %llx\n", cvt(RotationalVelocity[0]),
                cvt(RotationalVelocity[1]), cvt(RotationalVelocity[2]));
        fprintf(fp, "%llx %llx\n", cvt(NuclearEnergyTotal), cvt(AbsorbedEnergyTotal));
        fprintf(fp, "%lld\n", FlagDivideFile);
        fprintf(fp, "%lld %lld %lld %lld\n", FlagDamping, FlagNuclear, FlagBinary, FlagPotential);
        fprintf(fp, " %llx %llx\n", cvt(CodeUnit::BlackHoleMass),
                cvt(CodeUnit::BlackHoleMassInThisUnit));

        PS::S32 nproc = PS::Comm::getNumberOfProc();
        for(PS::S32 i = 0; i < nproc; i++) {
            PS::F64ort domain = dinfo.getPosDomain(i);
            fprintf(fp, "%llx %llx %llx %llx %llx %llx\n",
                    cvt(domain.low_[0]),  cvt(domain.low_[1]),  cvt(domain.low_[2]),
                    cvt(domain.high_[0]), cvt(domain.high_[1]), cvt(domain.high_[2]));
        }
    }

    void outputRunParameter(char **argv) {
        if(PS::Comm::getRank() != 0) {
            return;
        }
        FILE * fp = FilePointerForLog;
        fprintf(fp, "# # of process %8d\n", PS::Comm::getNumberOfProc());
        fprintf(fp, "# # of thread  %8d\n", PS::Comm::getNumberOfThread());
        if(atoi(argv[1]) == 0) {
            fprintf(fp, "# Header Data:  %s\n", argv[3]);
        } else {
            fprintf(fp, "# Restart file: %s\n", argv[2]);
        }
        fprintf(fp, "# CoulombFrc: %+e\n", CodeUnit::FractionOfCoulombCorrection);
        fprintf(fp, "# AlphaMax:   %+e\n", AlphaMaximum);
        fprintf(fp, "# AlphaMin:   %+e\n", AlphaMinimum);
        fprintf(fp, "# AlphuMax:   %+e\n", AlphuMaximum);
        fprintf(fp, "# AlphuMin:   %+e\n", AlphuMinimum);
        fprintf(fp, "# C for Time: %+e\n", CoefficientOfTimestep);
        fprintf(fp, "# KernelMax:  %+e\n", KernelSupportRadiusMaximum * CodeUnit::UnitOfLength);
        fprintf(fp, "# EpsilonU:   %+e\n", EpsilonOfInternalEnergy * CodeUnit::UnitOfEnergy);
        fprintf(fp, "# Damping: %d\n", FlagDamping);
        fprintf(fp, "# Nuclear: %d\n", FlagNuclear);        
        fprintf(fp, "# Binary:  %d\n", FlagBinary);        
        if(FlagBinary == 2) {
            fprintf(fp, "# BHPotential: %d\n", FlagPotential);
        }
        fprintf(fp, "# TimestepAscii: %+e\n", TimestepAscii);
        fprintf(fp, "# TimestepHexa:  %+e\n", TimestepHexa);
        fprintf(fp, "# # of Step %8d\n", NumberOfStep);
        fflush(fp);
    }

    template <class Tsph>
    PS::F64 setEpsilonOfInternalEnergy(Tsph & sph) {
        PS::F64 uminloc = std::numeric_limits<double>::max();
        for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            if(sph[i].uene < uminloc) {
                uminloc = sph[i].uene;
            }
        }
        PS::F64 uminglb = PS::Comm::getMinValue(uminloc);
        return 1e-4 * uminglb;
    }
    
}

template <class Tdinfo,
          class Tsph,
          class Tbhns,
          class Tmsls>
void readHexa(char * const filename,
              Tdinfo & dinfo,
              Tsph & sph,
              Tbhns & bhns,
              Tmsls & msls) {
    FILE *fp = fopen(filename, "r");
    RP::readHexa(fp, dinfo, sph, bhns, msls);
    for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].readHexa(fp);
    }
    for(PS::S32 i = 0; i < msls.getNumberOfParticleLocal(); i++) {
        msls[i].readHexa(fp);
    }
    for(PS::S32 i = 0; i < bhns.getNumberOfParticleLocal(); i++) {
        bhns[i].readHexa(fp);
    }
    fclose(fp);
    RP::NumberOfParticle = PS::Comm::getSum(sph.getNumberOfParticleLocal());
}

template <class Tdinfo,
          class Tsph,
          class Tbhns,
          class Tmsls>
void writeHexa(char * const filename,
               Tdinfo & dinfo,
               Tsph & sph,
               Tbhns & bhns,
               Tmsls & msls) {

    FILE *fp = fopen(filename, "w");
    RP::writeHexa(fp, dinfo, sph, bhns, msls);
    for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].writeHexa(fp);
    }
    for(PS::S32 i = 0; i < msls.getNumberOfParticleLocal(); i++) {
        msls[i].writeHexa(fp);
    }
    for(PS::S32 i = 0; i < bhns.getNumberOfParticleLocal(); i++) {
        bhns[i].writeHexa(fp);
    }
    fclose(fp);
}

template <class Tdinfo,
          class Tsph,
          class Tbhns,
          class Tmsls>
void outputData(Tdinfo & dinfo,
                Tsph & sph,
                Tbhns & bhns,
                Tmsls & msls) {
    if(RP::Time - (PS::S64)(RP::Time / RP::TimestepAscii) * RP::TimestepAscii == 0.) {
        char filename[64];
        if(RP::Time == 0.) {
            RP::NumberOfAscii = 0;
        }
        if(RP::FlagDivideFile == 0) {
            sprintf(filename, "snap/sph_t%04d.dat", RP::NumberOfAscii);
            sph.writeParticleAscii(filename);
        } else {
            sprintf(filename, "snap/sph_t%04d", RP::NumberOfAscii);
            sph.writeParticleAscii(filename, "%s_p%06d_i%06d.dat");
        }
        RP::NumberOfAscii++;
    }
    if(RP::Time - (PS::S64)(RP::Time / RP::TimestepHexa) * RP::TimestepHexa == 0.) {
        char filename[64];
        if(RP::Time == 0.) {
            RP::NumberOfHexa = 0;
        }
        sprintf(filename, "snap/t%04d_p%06d.hexa", RP::NumberOfHexa, PS::Comm::getRank());
        RP::NumberOfHexa++;
        writeHexa(filename, dinfo, sph, bhns, msls);
    }
    PS::F64 etot = calcEnergy(sph, bhns);
    PS::F64 enuc = 0.;
    PS::F64 vdet = calcDetonationVelocity(sph);
    if(PS::Comm::getRank() == 0) {
        using namespace CodeUnit;
        PS::F64 eabs = RP::AbsorbedEnergyTotal;
        fprintf(RP::FilePointerForLog, "time: %8d %16.10f %+e %+e %+e %+e %+e %+e %+e\n",
                RP::NumberOfStep, RP::Time * UnitOfTime, RP::Timestep * UnitOfTime,
                (etot - enuc - eabs) * UnitOfEnergy * UnitOfMass, etot * UnitOfEnergy * UnitOfMass,
                enuc * UnitOfEnergy * UnitOfMass, eabs * UnitOfEnergy * UnitOfMass,
                RP::KernelSupportRadiusMaximum * UnitOfLength, WT::getTimeTotal());
        fflush(RP::FilePointerForLog);
        WT::dump(RP::Time, RP::FilePointerForTime);
        fflush(RP::FilePointerForTime);
    }
    calcQuadrupoleMomentDot2(sph, bhns);
}

template <class Tsph,
          class Tbhns>
void predict(Tsph & sph,
             Tbhns & bhns) {
#ifdef USE_XEONPHI
#pragma omp parallel for
#endif
    for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].predict(RP::Timestep);
    }
    for(PS::S64 i = 0; i < bhns.getNumberOfParticleLocal(); i++) {
        bhns[i].predict(RP::Timestep);
    }
}

template <class Tsph,
          class Tbhns>
void correct(Tsph & sph,
             Tbhns & bhns) {
    if(RP::FlagDamping == 0) {
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            sph[i].correct(RP::Timestep);
        }
    } else if(RP::FlagDamping == 1) {
        for(PS::S64 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            sph[i].correctDamping1(RP::Timestep);
        }
    } else {
        fprintf(stderr, "Not supported damping mode %d!\n", RP::FlagDamping);
        PS::Abort();
    }
}

void initializeSimulation() {
    using namespace RunParameter;
    //Time              = 0.;
    //TimeEnd           = 0.;
    MaximumTimestep   = 10.;
    MinimumTimestep   = 1e-16;
    Timestep          = RP::MaximumTimestep;
    //TimeAscii         = 0.;
    //TimestepAscii     = 0.;
    //TimeHexa          = 0.;
    //TimestepHexa      = 0.;
    NumberOfStep      = 0;
    //NumberOfAscii     = 0;
    //NumberOfHexa      = 0;
    ComputationalBox.high_ = 0;
    ComputationalBox.low_  = 0;
    NumberOfDimension = 3;
    //KernelType        = 1;
    //CoefficientOfTimestep      = 0.1;
    GravitationalSoftening     = 3e6 * CodeUnit::UnitOfLengthInv;
    KernelSupportRadiusMaximum = 0.;
    EpsilonOfInternalEnergy    = 0.;    
    //NumberOfParticle   = 0;
    //AdiabaticIndex;
    RotationalVelocity = 0.;
    NuclearEnergyTotal = 0.; // enuc
    AbsorbedEnergyTotal = 0.;
    FlagGravity = 1;
    //FlagDamping = 0;
    //FlagNuclear = 0;
    if(PS::Comm::getRank() == 0) {
        RP::FilePointerForLog  = fopen("snap/time.log", "w");
        RP::FilePointerForTime = fopen("snap/prof.log", "w");
        RP::FilePointerForQuad = fopen("snap/quad.log", "w");
    }
    char debugfile[64];
    sprintf(debugfile, "snap/debug_p%06d.log", PS::Comm::getRank());
    RP::FilePointerForDebug = fopen(debugfile, "w");
}

template <class Tdinfo,
          class Tsph,
          class Tbhns,
          class Tmsls,
          class Tdensity,
          class Thydro,
          class Tgravity>
void startSimulation(char **argv,
                     Tdinfo & dinfo,
                     Tsph & sph,
                     Tbhns & bhns,
                     Tmsls & msls,
                     Tdensity & density,
                     Thydro & hydro,
                     Tgravity & gravity) {
    char headname[256];
    sprintf(headname, "%s.head", argv[3]);
    RP::readAscii(headname);
    if(atoi(argv[2]) == 0) {
        char filename[256];
        sprintf(filename, "%s.data", argv[3]);
        sph.readParticleAscii(filename);    
        RP::FlagDivideFile = 0;
    } else {
        sph.readParticleAscii(argv[3], "%s_p%06d_i%06d.data");    
        RP::FlagDivideFile = 1;
    }
    {
        RP::FlagBinary = 0;
        char filename[256];
        sprintf(filename, "%s.bhns", argv[3]);
        FILE *fp  = fopen(filename, "r");
        if(fp != NULL) {
            fclose(fp);
            bhns.readParticleAscii(filename);
            RP::FlagBinary = 1;
        }
        sprintf(filename, "%s.imbh", argv[3]);
        fp = fopen(filename, "r");
        if(fp != NULL) {
            using namespace CodeUnit;
            RP::FlagBinary = 2;
            fscanf(fp, "%lf%lld", &BlackHoleMass, &RP::FlagPotential);
            fclose(fp);
            BlackHoleMassInThisUnit = BlackHoleMass * SolarMass * UnitOfMassInv;
        }
    }
    ND::setDimension(RP::NumberOfDimension);
    SK::setKernel(RP::KernelType, RP::NumberOfDimension);
    RP::KernelSupportRadiusMaximum = std::numeric_limits<PS::F64>::max();
    RP::EpsilonOfInternalEnergy    = RP::setEpsilonOfInternalEnergy(sph);
    RP::outputRunParameter(argv);
    PS::MT::init_genrand(0);
    dinfo.decomposeDomainAll(sph);
    sph.exchangeParticle(dinfo);
    bhns.exchangeParticle(dinfo);
    calcSPHKernel(dinfo, sph, bhns, msls, density, hydro, gravity);
}

template <class Tdinfo,
          class Tsph,
          class Tbhns,
          class Tmsls>
void restartSimulation(char **argv,
                       Tdinfo & dinfo,
                       Tsph & sph,
                       Tbhns & bhns,
                       Tmsls & msls) {
    char filename[64];
    sprintf(filename, "%s_p%06d.hexa", argv[2], PS::Comm::getRank());
    readHexa(filename, dinfo, sph, bhns, msls);
    ND::setDimension(RP::NumberOfDimension);
    SK::setKernel(RP::KernelType, RP::NumberOfDimension);
    // *********************
    //RP::CoefficientOfTimestep = 1e-2;
    //RP::CoefficientOfTimestep = 1e-1;
    //RP::TimestepAscii   = 1. / 256.;
    //RP::MaximumTimestep = 1. / 1024.;
    //RP::KernelSupportRadiusMaximum *= 0.5;
    //RP::ReductionTime *= 2.0;
    // *********************
    RP::outputRunParameter(argv);
    return;
}

template <class Tsph>
void initializeWriteTemperatureDensity(Tsph & sph) {
    FILE *fp = fopen("id.list", "r");
    if(fp != NULL) {
        PS::S64 nlist = 0;
        fscanf(fp, "%lld", &nlist);
        PS::S64 *list = (PS::S64 *)malloc(sizeof(PS::S64) * nlist);
        for(PS::S64 j = 0; j < nlist; j++) {
            fscanf(fp, "%lld", &list[j]);
        }
        for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
            for(PS::S64 j = 0; j < nlist; j++) {
                if(sph[i].id == list[j]) {
                    sph[i].flagwrite = true;
                    break;
                }
            }
        }
        free(list);
    }
    return;
}

template <class Tsph>
void writeTemperatureDensity(Tsph & sph) {
    for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        sph[i].writeTemperatureDensity(RP::FilePointerForDebug);
    }
    fflush(RP::FilePointerForDebug);
}

template <class Tsph>
void dumpOneParticle(PS::S32 id,
                     Tsph & sph) {
    for(PS::S32 i = 0; i < sph.getNumberOfParticleLocal(); i++) {
        if(sph[i].id == id) {
            fprintf(RP::FilePointerForDebug, "%16.10f ", RP::Time * CodeUnit::UnitOfTime);
            sph[i].writeAscii(RP::FilePointerForDebug);
        }
    }
    fflush(RP::FilePointerForDebug);
}

template <class Tdinfo,
          class Tsph,
          class Tbhns,
          class Tmsls,
          class Tdensity,
          class Thydro,
          class Tgravity>
void loopSimulation(Tdinfo & dinfo,
                    Tsph & sph,
                    Tbhns & bhns,
                    Tmsls & msls,
                    Tdensity & density,
                    Thydro & hydro,
                    Tgravity & gravity) {
    initializeWriteTemperatureDensity(sph);
    WT::clear();
    while(RP::Time < RP::TimeEnd) {
        WT::reduceInterProcess();
        outputData(dinfo, sph, bhns, msls);
        WT::clear();
        RP::Timestep = calcTimestep(sph);
        WT::start();        
        WT::accumulateCalcNuclearReaction();
        WT::start();
        predict(sph, bhns);
        WT::accumulateOthers();
        WT::start();
        if(RP::NumberOfStep % 4 == 0) {
            PS::MT::init_genrand(0);
            dinfo.decomposeDomainAll(sph);
        }
        WT::accumulateDecomposeDomain();
        WT::start();
        sph.exchangeParticle(dinfo);
        bhns.exchangeParticle(dinfo);
        WT::accumulateExchangeParticle();
        calcSPHKernel(dinfo, sph, bhns, msls, density, hydro, gravity);
        WT::start();
        correct(sph, bhns);
        WT::accumulateOthers();
        RP::Time += RP::Timestep;
        RP::NumberOfStep++;
    }
}

template <class Tdinfo,
          class Tsph,
          class Tbhns,
          class Tmsls>
void finalizeSimulation(Tdinfo & dinfo,
                        Tsph & sph,
                        Tbhns & bhns,
                        Tmsls & msls) {
    outputData(dinfo, sph, bhns, msls);
    if(PS::Comm::getRank() == 0) {
        fclose(RP::FilePointerForLog);
        fclose(RP::FilePointerForTime);
        fclose(RP::FilePointerForQuad);
    }
    fclose(RP::FilePointerForDebug);
    PS::F64    mc;
    PS::F64vec xc;
    PS::F64vec vc;
    if(RP::FlagBinary == 0) {
        calcCenterOfMass(sph, mc, xc, vc);
    } else if(RP::FlagBinary == 1) {
        calcCenterOfMass(sph, bhns, mc, xc, vc);
    } else {
        ;
    }
    PS::S32 nloc = sph.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nloc; i++) {
        sph[i].pos -= xc;
        sph[i].vel -= vc;
    }
    if(RP::FlagDamping == 1) {
	for(PS::S32 i = 0; i < nloc; i++) {
	    PS::F64vec vrot = RP::RotationalVelocity ^ sph[i].pos;
	    sph[i].vel += vrot;
	}
    }
    if(RP::FlagDivideFile == 0) {
        sph.writeParticleAscii("snap/final.dat");
    } else {
        char filename[64];
        sprintf(filename, "snap/final");
        sph.writeParticleAscii(filename, "%s_p%06d_i%06d.dat");
    }
}

typedef HelmholtzGas GeneralSPH;

