#pragma once



namespace nhfInt {

class EriController {
public:
    bool    eriScreenSwitch;
    double  eriScreenCutoff;

    EriController():
    eriScreenSwitch(true),
    eriScreenCutoff(1e-8)  {}
};

} // namespace (nhfInt)