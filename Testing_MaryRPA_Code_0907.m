NumRepetitions = 80;

RPA = 30;
LengtH = 1000;
DNA = zeros(1,LengtH);


% initial saturation of ssDNA with RPA %
for n = 1:NumRepetitions 

        % chooses point of attachment and gives random probabilities for RPA binding %
        RRPA = randi((LengtH-RPA),1);
        RPAProbability = rand;
        
            % RPA binding: first half of RPA (labeled A) = 4, second half (labeled D) = 16 %
            if (DNA(RRPA:RRPA+(RPA-1))) == 0 & RPAProbability <= 0.5 %#ok<AND2>
            DNA(RRPA:RRPA+((RPA/2)-1)) = 4;
            DNA((RRPA+(RPA/2)):(RRPA+(RPA-1))) = 16;
            end
end

RPA_Count = (length(find(DNA == 4))+length(find(DNA == 16)))/RPA;
Saturation = (RPA_Count*RPA)/LengtH;