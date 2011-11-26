#ifndef ROOT_ToyMCSamplerOpt_h
#define ROOT_ToyMCSamplerOpt_h

#include <RooStats/ToyMCSampler.h>
struct RooProdPdf;
struct RooPoisson;

namespace toymcoptutils {
    class SinglePdfGenInfo {
        public:
            enum Mode { Binned, Poisson, Unbinned, Counting };
            SinglePdfGenInfo(RooAbsPdf &pdf, const RooArgSet& observables, bool preferBinned, const RooDataSet* protoData = NULL, int forceEvents = 0) ;
            ~SinglePdfGenInfo() ;
            RooAbsData *generate(const RooDataSet* protoData = NULL, int forceEvents = 0) ;
            RooDataSet *generateAsimov(RooRealVar *&weightVar) ;
            const RooAbsPdf * pdf() const { return pdf_; }
        private:
            Mode mode_;
            RooAbsPdf *pdf_; 
            RooArgSet observables_;
            RooAbsPdf::GenSpec *spec_;
            RooRealVar *weightVar_;
            RooDataSet *generateWithHisto(RooRealVar *&weightVar, bool asimov) ;
            RooDataSet *generateCountingAsimov() ;
            void setToExpected(RooProdPdf &prod, RooArgSet &obs) ;
            void setToExpected(RooPoisson &pois, RooArgSet &obs) ;
    };
    class SimPdfGenInfo {
        public:
            SimPdfGenInfo(RooAbsPdf &pdf, const RooArgSet& observables, bool preferBinned, const RooDataSet* protoData = NULL, int forceEvents = 0) ;
            ~SimPdfGenInfo() ;
            RooAbsData *generate(RooRealVar *&weightVar, const RooDataSet* protoData = NULL, int forceEvents = 0) ;
            RooAbsData *generateAsimov(RooRealVar *&weightVar) ;
            void setCopyData(bool copyData) { copyData_ = copyData; }
        private:
            RooAbsPdf                       *pdf_; 
            RooAbsCategoryLValue            *cat_;
            RooArgSet                        observables_;
            std::vector<SinglePdfGenInfo *>  pdfs_; 
            RooArgSet                        ownedCrap_;
            std::map<std::string,RooAbsData*> datasetPieces_;
            bool                              copyData_;
            //std::map<std::string,RooDataSet*> datasetPieces_;

    }; 
}

class ToyMCSamplerOpt : public RooStats::ToyMCSampler{
    public:
        ToyMCSamplerOpt(RooStats::TestStatistic& ts, Int_t ntoys, RooAbsPdf *globalObsPdf = 0) ;
        ToyMCSamplerOpt(const RooStats::ToyMCSampler &base) ;
        ToyMCSamplerOpt(const ToyMCSamplerOpt &other) ;
        ~ToyMCSamplerOpt() ;
        virtual RooAbsData* Generate(RooAbsPdf& pdf, RooArgSet& observables, const RooDataSet* protoData = NULL, int forceEvents = 0) const ;
        virtual void SetPdf(RooAbsPdf& pdf) ;
        void setGlobalObsPdf(RooAbsPdf *pdf) { globalObsPdf_ = pdf; }
    private:
        RooAbsPdf *globalObsPdf_;
        mutable RooDataSet *globalObsValues_; 
        mutable int globalObsIndex_;

        mutable RooRealVar *weightVar_;
        mutable std::map<RooAbsPdf *, toymcoptutils::SimPdfGenInfo *> genCache_;
#if ROOT_VERSION_CODE < ROOT_VERSION(5,29,0)
    public:
        virtual RooAbsData* GenerateToyData(RooArgSet& /*nullPOI*/) const ;
    private:
        // objects below cache information and are mutable and non-persistent
        mutable RooArgSet* _allVars ; //! 
        //mutable std::list<RooAbsPdf*> _pdfList ; //!
        //mutable std::list<RooArgSet*> _obsList ; //!
        //mutable std::list<RooAbsPdf::GenSpec*> _gsList ; //!      
#endif
        
};

#endif
