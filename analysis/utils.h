#ifndef UTILS
#define UTILS 1

void ScaleHistogramsToGlobalMax(TH2* histograms[], int size) {
    if (size <= 0) {
        std::cerr << "Error: The array size must be greater than zero!" << std::endl;
        return;
    }

    double globalMax = -std::numeric_limits<double>::max();
    double globalMin =  std::numeric_limits<double>::max();

    // Step 1: Find the true global min/max from bin contents
    for (int i = 0; i < size; ++i) {
        if (!histograms[i]) continue;

        // Reset any user-defined scale so Get{Max,Min}imum()
        // recompute from the actual bin contents (key fix after Clone/Add)
        histograms[i]->SetMaximum();   // -> -1111 ("unset")
        histograms[i]->SetMinimum();   // -> -1111 ("unset")

        double localMax = histograms[i]->GetMaximum();
        double localMin = histograms[i]->GetMinimum();

        if (localMax > globalMax) globalMax = localMax;
        if (localMin < globalMin) globalMin = localMin;   // fixed: '<' not '>'
    }

    if (globalMax <= globalMin) {
        std::cerr << "Error: invalid global range (max <= min)!" << std::endl;
        return;
    }

    // Step 2: Apply the common scale to every histogram
    for (int i = 0; i < size; ++i) {
        if (histograms[i]) {
            histograms[i]->SetMaximum(globalMax);
            histograms[i]->SetMinimum(globalMin);
        }
    }

    std::cout << "Histograms scaled to global range ["
              << globalMin << ", " << globalMax << "]" << std::endl;
}

void SetDivergingPalette(int nContours = 255) {
    const int nStops = 3;

    // position of each color stop along [0,1]
    Double_t stops[nStops] = { 0.00, 0.50, 1.00 };

    //                         neg     zero    pos
    Double_t red[nStops]   = { 0.13,  1.00,  0.70 };
    Double_t green[nStops] = { 0.27,  1.00,  0.05 };
    Double_t blue[nStops]  = { 0.65,  1.00,  0.15 };

    TColor::CreateGradientColorTable(nStops, stops, red, green, blue, nContours);
    gStyle->SetNumberContours(nContours);
}

#endif
