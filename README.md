# A Quadratic Programming Algorithm with $O(n^3)$ Time Complexity 
This is the code repository for https://arxiv.org/pdf/2507.04515

Our manuscript has four main contributions:

(1) demonstrates how to transform a general strictly convex QP, a Lasso problem, and a support vector machine (or regression) into a Box-QP, highlighting the broad applicability of the approach.

(2) develops a unified and simple proof framework for feasible IPM algorithms with \textit{exact Newton step} and \textit{approximated Newton step}, respectively, resulting in $O(n^{3.5})$ and $O(n^3)$ time complexity.

(3)  proves that the proposed $O(n^3)$-time-complexity algorithm has an \textit{exact} and \textit{data-independent} number of iterations   

<img width="273" alt="image" src="https://github.com/user-attachments/assets/13198375-3a04-4b80-a09f-74d6a28faadb" />

and the \textit{data-independent} total number of rank-1 updates bounded by

<img width="249" alt="image" src="https://github.com/user-attachments/assets/3e6f768b-e4e8-4323-9335-2954d0009aa4" />

where $\alpha,\beta,\eta,\delta$ are \textit{data-independent} constants, thereby being able to offer an execution time certificate for real-time optimization-based applications;

(4) shows that the proposed $O(n^3)$-time-complexity algorithm is simple to implement and matrix-free in the iterative procedures.

## Numerical Validations
#### (1) Validation of number of iterations and rank-1 updates
<img width="479" alt="image" src="https://github.com/user-attachments/assets/41505dab-a47d-45cc-b2e4-5b76d2100779" />

<img width="487" alt="image" src="https://github.com/user-attachments/assets/395e9bce-d0d3-4716-b26f-1e24eb76fbb7" />

#### (2) Numerical comparison of Algorithms 1 ($O(n^{3.5)}$) and 2 ($O(n^{3)}$)
<img width="495" alt="image" src="https://github.com/user-attachments/assets/6bfd7c0a-39ef-4ab5-8afb-bc4ba5b1a922" />


