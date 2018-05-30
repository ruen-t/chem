angular
    .module('app')
    .controller('DialogController', DialogController);
function DialogController($scope, $mdDialog, $http) {
    var self = this;
    self.description = "";
    self.create = create;
    $scope.$watch("vm.smart", function (newVal, oldVal) {
        self.isCorrectSmart = true;
    });
    $scope.hide = function () {
        $mdDialog.hide();
    };

    $scope.cancel = function () {
        $mdDialog.cancel();
    };

    function create() {
        isCorrectSmart().then(r => {
            if (self.isCorrectSmart) {
                $mdDialog.hide();
                var name = self.name;
                var smart = self.smarts;
                var description = self.description;
                var payload = new FormData();
                var reaction = { name, smart, description };
                payload.append("reaction", JSON.stringify(reaction));
                $http({
                    url: REACTION_API,
                    method: 'POST',
                    data: payload,
                    headers: { 'Content-Type': undefined },
                    transformRequest: angular.identity
                }).then(function (response) {
                    console.log(response);
                    swal('Successfully add new reaction', '', 'success')
                }, function (error) {
                    swal('Fail to add new reaction', '', 'error')
                })
            }
        })
    }
    function isCorrectSmart() {
        return new Promise(resolve => {
            $http({
                url: CHECK_SMART_API + self.smarts,
                method: 'GET',
                headers: { 'Content-Type': 'application/json' },
                transformRequest: angular.identity
            }).then(function (response) {
                console.log(response);
                self.isCorrectSmart = response.data.result;
                resolve(response.data.result);
            }, function (error) {
                console.log(error)
                self.isCorrectSmart = false;
                resolve(false);
            })
        })

    }
}