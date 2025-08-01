# TITE-STEIN: A Time-to-Event Model-Assisted Design for Phase I/II Trials

**TITE-STEIN** (Time-to-event Simple Toxicity and Efficacy Interval Design) is a model-assisted dose-finding design developed to accelerate Phase I/II clinical trials by incorporating time-to-event toxicity and efficacy outcomes. It extends the STEIN framework to handle late-onset outcomes and introduces an **Optimal Biological Dose (OBD) verification** procedure to mitigate the risk of selecting inadmissible doses.

## Features
- Supports real-time dose decision-making using partial follow-up information
- Incorporates both toxicity and efficacy outcomes for OBD selection
- Includes a verification step to avoid selecting suboptimal or non-existent OBDs
- Demonstrated improved operating characteristics compared to other TITE-based designs through extensive simulations
- Straightforward to implement with no complex dose-toxicity and dose-efficacy modeling assumptions

## Repository Contents

- `.R`: Main function used for simulation studies in *Sun, H., Tu, J., Ananthakrishnan, R., & Kim, E. (2025). TITE-STEIN: Time-to-event Simple Toxicity and Efficacy Interval Design to Accelerate Phase I/II Trials. Journal of Biopharmaceutical Statistics (under review).*
- `.R`: 

## Getting Started
```
# Load source files

# Define parameters (e.g., true toxicity and efficacy probabilities for 5 dose levels)

# Run the trial simulation

# View OBD selection result

```
## Reference
If you use this code, please cite: 
> Sun, H., Tu, J., Ananthakrishnan, R., & Kim, E. (2025). TITE-STEIN: Time-to-event Simple Toxicity and Efficacy Interval Design to Accelerate Phase I/II Trials. *Journal of Biopharmaceutical Statistics* (under review).
