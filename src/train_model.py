import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score
from sklearn.preprocessing import LabelEncoder

# 1) Load dataset
df = pd.read_csv("D:/project_fastq/ml_dataset.csv")

print("Dataset shape:", df.shape)

# 2) Encode labels (stroke / healthy)
le = LabelEncoder()
df["disease"] = le.fit_transform(df["disease"])

# stroke = 1
# healthy = 0

# 3) Split features and label
X = df.drop(columns=["disease", "sample"], errors="ignore")
y = df["disease"]

# 4) Train/test split
X_train, X_test, y_train, y_test = train_test_split(
    X, y,
    test_size=0.2,
    random_state=42,
    stratify=y
)

# 5) Model
model = RandomForestClassifier(
    n_estimators=200,
    random_state=42
)

# 6) Train
model.fit(X_train, y_train)

# 7) Predict
y_pred = model.predict(X_test)

# 8) Evaluation
print("\nAccuracy:", accuracy_score(y_test, y_pred))
print("\nClassification Report:\n")
print(classification_report(y_test, y_pred, target_names=le.classes_))