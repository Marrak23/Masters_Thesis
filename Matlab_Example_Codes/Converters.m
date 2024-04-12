
%--------------------Converters--------------%
classdef Converters
    methods
        %---------Energy.---------------------------%
        %From eV to J---------------%
        function ToJ = FromeVtoJ(~, FromeV)
            ToJ = FromeV * 1.6022*10^(-19);
        end
        %---------------------------------------%

        %From J to eV---------------------------%
        function ToeV = FromJtoeV(~, FromJ)
            ToeV = FromJ*6.242*10^(18);
        end
        %---------------------------------------%
        
        %From eV to Hartree---------------------------%
        function ToHartree = FromeVtoHartree(~, FromeV)
            ToHartree =  FromeV*0.03674930495120813;
        end
        %---------------------------------------%

        %From Hartree to eV---------------------------%
        function ToeV = FromHartreetoeV(~, FromHartree)
            ToeV =  FromHartree*(0.03674930495120813)^(-1);
        end
        %---------------------------------------%

        %---------Energy.---------------------------%


        %From A to M---------------------------%
        function ToM = FromAtoM(~, FromA)
            ToM = FromA*10^(-10);
        end
        %---------------------------------------%


        %From kg to electron kg-----------------%
        function ToEkg = FromKgtoEkg(~, FromKg)
            ToEkg =FromKg*(9.1093837*10^(-31))^(-1);
        end
        %---------------------------------------%

        %From kg to electron kg-----------------%
        function Tokg = FromEkgtoKg(~, FromEkg)
            Tokg =FromEkg*9.1093837*10^(-31);
        end
        %---------------------------------------%
        
        %From A to a0-----------------%
        function Toa0 = FromAtoa0(~, FromA)
            Toa0 = FromA*1.88973;
        end
        %---------------------------------------%

        %From a0 to A-----------------%
        function ToA = Froma0toA(~, Froma0)
            ToA = Froma0*(1.88973)^(-1);
        end
        %---------------------------------------%
        
        %From M to A----------------------------%
        function ToA = FromMtoA(~, FromM)
            ToA = FromM*10^(10);
        end
        %---------------------------------------%
            
        %From m to cm---------------------------%
        function ToCm = FromMtoCm(~, FromM)
            ToCm = FromM*10^(2);
        end
        %---------------------------------------%

        %From cm to m---------------------------%
        function ToM = FromCmtoM(~, FromCm)
            ToM = FromCm*10^(-2);
        end
        %---------------------------------------%

        %From hbar/Ehartree to femtosecond------&
            function ToFs = FromhbarEharttoFs(~, FromhbarEhart)
                ToFs = FromhbarEhart*0.02418884254;
            end
        %---------------------------------------%
        
        %From femtosecond to hbar/Ehartree------%
            function TohbarEhart = FromFstohbarEhart(~, FromFs)
                TohbarEhart = FromFs*0.02418884254^(-1);
            end
        %---------------------------------------%
        
    end
end
%---------------------------------------%




