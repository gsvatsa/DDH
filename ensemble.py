from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import VotingRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error
import matplotlib.pyplot as plt

def get_ensemble_model2(X, y):
  X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
  reg1 = GradientBoostingRegressor(n_estimators=100, learning_rate=0.1, max_depth=1, random_state=1, loss='ls').fit(X_train, y_train)
  mean_squared_error(y_test, reg1.predict(X_test))
  reg2 = RandomForestRegressor(random_state=1, n_estimators=100).fit(X_train, y_train)
  mean_squared_error(y_test, reg2.predict(X_test))
  reg3 = LinearRegression().fit(X_train,y_train)
  mean_squared_error(y_test, reg3.predict(X_test))
  ereg = VotingRegressor(estimators=[('gb', reg1), ('rf', reg2), ('lr', reg3)]).fit(X_train, y_train)
  mean_squared_error(y_test, ereg.predict(X_test))
  return ereg

# X = r.y_hats.iloc[:, 0:2]
# y = r.y_hats.iloc[:,2]
# ereg = get_ensemble_model2(X, y)
# mean_squared_error(y, ereg.predict(X))
# pred1 = reg1.predict(X)
# pred2 = reg2.predict(X)
# pred3 = reg3.predict(X)
# pred4 = ereg.predict(X)
# 
# plt.figure()
# plt.plot(pred1, 'gd', label='GradientBoostingRegressor')
# plt.plot(pred2, 'b^', label='RandomForestRegressor')
# plt.plot(pred3, 'ys', label='LinearRegression')
# plt.plot(pred4, 'r*', ms=10, label='VotingRegressor')
# 
# plt.tick_params(axis='x', which='both', bottom=False, top=False,
#                 labelbottom=False)
# plt.ylabel('predicted')
# plt.xlabel('training samples')
# plt.legend(loc="best")
# plt.title('Regressor predictions and their average')
# 
# plt.show()



