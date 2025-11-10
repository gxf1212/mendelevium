

4.1.1 Regression Metrics
For regression models, there are several metrics (Table 6.4) that should be calculated
to assess both robustness and predictability. Robustness could be defined as the
ability to maintain its predictive performance across different data sets and is
Fig. 6.6 Validations and controls for training and test/repeated k-Fold/Leave-One-Out (LOO)
144 P. O. Fernandes and V. G. Maltarollocalculated using different internal validation strategies. The predictability is the
predictive performance of a QSAR model for unseen new molecules (test set). For
example, R2-derived metrics, such as Q2 calculated with LOO or other internal
validation strategy, a metric to estimate the robustness of a model. Usually, some
authors report an R2 (calibration coefficient) metric which is calculated with the
same formula as Q2 but with no internal validation. Together, R2 and Q2 could be
used as a guide to prevent overfitting: ideally, the difference between those two
metrics should be lower than 0.3 (R2 – Q2 > 0.3) [93, 94].
In addition, errors such as RMSE and/or MAE metrics could be calculated in both
internal and external validation procedures. The Root Mean Squared Error (RMSE),
derived from MSE, offers an interpretable measure by maintaining the same units as
the target variable. The Mean Absolute Error (MAE), which measures the average
magnitude of errors in predictions without considering their direction, provides a
straightforward interpretation of prediction accuracy.
Nowadays, it is recommended that a consensus of all metrics should be used in
the evaluation of a given QSAR model’s quality [110]. Those metrics are Q2
coefficients discussed by Consonni and colleagues [112], r2m parameters were introduced by Roy and collaborators [111], and CCC from Gramatica work [113] as well
as the errors of predictions.

